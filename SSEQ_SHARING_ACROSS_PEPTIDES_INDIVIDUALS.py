"""
scripts to parse SSEQ TCR data from HIMC using Python.
@author: Aditya Ambati, Mignot Lab, Stanford University
"""
import sys
import argparse
import datetime
import os
import pandas as pd 
from collections import defaultdict

parser = argparse.ArgumentParser(description='A class method to process sseq libraries to pick public cdr3s')
parser.add_argument('-file', required=True, help='Full path to sseq library files')
args=parser.parse_args()
filein = args.file

####### define the sseq as an object ########
class sseq_calls(object):
	"""A sseq seqeuncing object:

	Attributes:
		name: A string representing the sseq library.
		filein: location of the associated sseq library in csv format.
	"""
	instance_count =0

	def __init__(self, name, filein):
		"""intitate the sseq object class."""
		self.name = name
		self.filein = filein
		self.instance_count += 1

	# Not used		
	def get_attr(self):
		print(' NAME :- {} \n SOURCE FILE :- {} \n INSTANCE COUNT :- {}'.format(self.name, self.filein, self.instance_count))

	
	def read_data(self):
		data=[]
		with open(self.filein) as sseq_in:
			print('READING DATA FROM FILE  {} '.format(self.filein))
			for line in sseq_in:
				data.append(line.strip('\n'))

		return data

	@staticmethod
	def norm_fam_cdrs(fam_dic):
		"""will take a cdr3 and associated families and normalize the fam names for same cdrs based on the counts"""
		norm_cdrs =defaultdict()
		for k, v in fam_dic.items():
			if len(v) > 1:
				norm_cdrs[k] = max(v.items(), key = lambda a: a[0])[0]
			else: 
				norm_cdrs[k] = list(v.keys())[0]

		return norm_cdrs

	#@classmethod
	def process_data_CDR3(self):

		# Initialize outputs
		fam_dic_alpha=defaultdict(lambda: defaultdict(int))
		fam_dic_beta=defaultdict(lambda: defaultdict(int))
		CDR3_DIC=defaultdict(lambda: defaultdict(int))
		CDR3a_DIC =defaultdict(lambda: defaultdict(int))
		CDR3b_DIC =defaultdict(lambda: defaultdict(int))
		CDR3_totals =defaultdict(lambda: defaultdict(int))
		data_list = self.read_data()

		# For every instance
		for n, line in enumerate(data_list):
			if ";" in line:
				parse_line=line.strip('\n').split(';')
			else:
				parse_line=line.strip('\n').split(',')
			if parse_line[1] and parse_line[2] and parse_line[3] and parse_line[0]:
				make_key = ','.join([parse_line[1], parse_line[2].strip(), parse_line[3]])

			VB=parse_line[6]
			VA=parse_line[12]
			JA = parse_line[14].split(' ')[0]
			VAalt = parse_line[18]
			JAalt = parse_line[19].split(' ')[0]
			
			cdr3a = parse_line[15]
			cdr3b = parse_line[9]
			cdr3a_alt= parse_line[20]

			if cdr3a and cdr3b and not cdr3a_alt:
				################## here this script picks up paired chains ##################################################
				#cdr=VA+'_'+cdr3a+'_'+JA+':'+VB+'_'+cdr3b
				cdr = cdr3a+':'+cdr3b
				CDR3_DIC[cdr][make_key] += 1
				CDR3_totals[make_key]['both'] += 1
			elif cdr3a_alt and cdr3b and not cdr3a:
				# cdr = VAalt+'_'+cdr3a_alt+'_'+JAalt+':'+VB+'_'+cdr3b
				cdr = cdr3a_alt+':'+cdr3b
				CDR3_DIC[cdr][make_key] += 1
				CDR3_totals[make_key]['both'] += 1
			elif cdr3a_alt and cdr3a and cdr3b:
				# cdr1 = VAalt+'_'+cdr3a_alt+'_'+JAalt+':'+VB+'_'+cdr3b
				# cdr2 = VA+'_'+cdr3a+'_'+JA+':'+VB+'_'+cdr3b
				cdr1 = cdr3a_alt+':'+cdr3b
				cdr2 =cdr3a+':'+cdr3b
				#print "FOUND alt alphas ALPHA {}, alt ALPHA {} BETAS {}".format(cdr3a, cdr3a_alt, cdr3b)
				for cdr in [cdr1, cdr2]:
					CDR3_DIC[cdr][make_key] += 1
					CDR3_totals[make_key]['both'] += 1

			if cdr3a and not cdr3a_alt:######################################################### here this script picks up only alpha chains #############
				#cdr3a = cdr3a
				CDR3a_DIC[cdr3a][make_key] += 1
				CDR3_totals[make_key]['a'] += 1
				fam_dic_alpha[cdr3a][VA+'_'+JA] += 1

			elif cdr3a_alt and not cdr3a:
				#cdr3a_alt= VAalt+'_'+cdr3a_alt+'_'+JAalt
				CDR3a_DIC[cdr3a_alt][make_key] += 1
				CDR3_totals[make_key]['a'] += 1
				fam_dic_alpha[cdr3a_alt][VAalt+'_'+JAalt] += 1
				
			elif cdr3a_alt and cdr3a:
				for cdr, V, J in zip([cdr3a_alt, cdr3a], [VAalt, VA], [JAalt, JA]):
					#cdr = V + '_' +cdr+'_'+J
					CDR3a_DIC[cdr][make_key] += 1
					CDR3_totals[make_key]['a'] += 1
					fam_dic_alpha[cdr][V+'_'+J] += 1

			if cdr3b:######################################################### here this script picks up only beta chains #############
				#cdr3b = VB+'_'+cdr3b
				CDR3b_DIC[cdr3b][make_key] += 1
				CDR3_totals[make_key]['b'] += 1
				fam_dic_beta[cdr3b][VB] += 1
	
		print('Processed {} single seq calls '.format(n))
		return [CDR3_DIC, CDR3a_DIC, CDR3b_DIC, CDR3_totals, fam_dic_beta, fam_dic_alpha]
	
	@staticmethod
	def get_normalized_fam(cdr, normfamdicalpha, normfamdicbeta, chain):
		"""fucntion returns a normalized familiy of a given cdr3"""
		#global make_k
		if chain == 'CDR3AB':
			get_alpha_CDR = cdr.split(':')[0]
			get_beta_CDR = cdr.split(':')[1]
			if get_alpha_CDR in normfamdicalpha:
				get_fam = normfamdicalpha[get_alpha_CDR].split('_')
				make_k1 = get_fam[0]+'_'+get_alpha_CDR+'_'+get_fam[1]
			if get_beta_CDR in normfamdicbeta:
				get_fam = normfamdicbeta[get_beta_CDR]
				make_k2 = get_fam+'_'+get_beta_CDR
			make_k = make_k1+':'+make_k2
		elif chain == 'CDR3A':
			if cdr in normfamdicalpha:
				get_fam = normfamdicalpha[cdr].split('_')
				make_k = get_fam[0]+'_'+cdr+'_'+get_fam[1]
		else:
			if cdr in normfamdicbeta:
				get_fam = normfamdicbeta[cdr]
				make_k = get_fam+'_'+cdr
				
		return make_k

	@staticmethod
	def keep_enriched_clones(CDR_DIC, depth, per_ind_depth):
		"""# takes a master dictionary and processses with a condition that a cdr3 must be present atleast in 3 individuals and more than 3x returns a keep list of cdrs"""
		out_cdr=defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
		for k, v in CDR_DIC.items():
			temp_dic ={}
			if len(v) > per_ind_depth:
				for k2, v2 in v.items():
					parse_k2 =  k2.split(',')
					dbid = parse_k2[0]+','+parse_k2[-1]
					pep = parse_k2[1]
					out_cdr[k][pep][dbid] += v2

		keep_list = []					
		for k, v in out_cdr.items():
			for k2, v2 in v.items():
				if len(v2) >= depth:
					#print k, k2
					if k not in keep_list:
						keep_list.append(k)
		return keep_list


	def write_outfile(self, cdr3dic, cdr3keep, normfamdicalpha, normfamdicbeta, chain, CDR3_totals):
		""" this function writes out the sharing results """
		
		dttime=datetime.datetime.now().strftime ("%Y%m%d")


		header=',DBID,PEPTIDE,DX,COUNTS,Totalcalls'+'\n'
		if chain == 'CDR3AB':
			outfile = open('CDR3AB_SHARING_SSEQ_'+dttime+'.csv', 'w')
			outfile.write('CDR3AB')
			cdr3_total_key = 'both'
		elif chain =='CDR3A':
			outfile = open('CDR3A_SHARING_SSEQ_'+dttime+'.csv', 'w')
			outfile.write('CDR3A')
			cdr3_total_key = 'a'
		else:
			outfile = open('CDR3B_SHARING_SSEQ_'+dttime+'.csv', 'w')
			outfile.write('CDR3B')
			cdr3_total_key = 'b'
		outfile.write(header)
		for k, v in cdr3dic.items():
			if k in cdr3keep:
				make_k=self.get_normalized_fam(cdr = k, normfamdicbeta=normfamdicbeta, normfamdicalpha=normfamdicalpha, chain=chain)
				if len(v) > 1:
					for k2, v2 in v.items():
						if k2 in CDR3_totals:
							if cdr3_total_key in CDR3_totals.get(k2):
								get_counts = CDR3_totals.get(k2).get(cdr3_total_key)
							if 'Blank,Blank,Blank' not in k2:
								outfile.write(make_k+','+k2+','+str(v2)+','+str(get_counts)+','+'\n')
		outfile.close()


	def conclude_process(self, depth_across_peptide, depth_per_individual):
		"""This will conclude the analysis writing out the sharing results"""

		processed_list = self.process_data_CDR3()
		if len(processed_list) == 6:
			assert len(processed_list) == 6
			CDR3_DIC, CDR3a_DIC, CDR3b_DIC, CDR3_totals, fam_dic_beta, fam_dic_alpha = [item for item in processed_list]
			print('ALL CDR3 SSEQ DATA HAS BEEN PROCESSED ')
			print('DEPTH OF THE SSEQ ACROSS PEPTIDE {} & WITHIN_INDIVIDUAL {}'.format(depth_across_peptide, depth_per_individual))
			cdr3a_keep=self.keep_enriched_clones(CDR_DIC = CDR3a_DIC, depth=depth_across_peptide, per_ind_depth=depth_per_individual)
			cdr3b_keep=self.keep_enriched_clones(CDR_DIC = CDR3b_DIC, depth=depth_across_peptide, per_ind_depth=depth_per_individual)
			cdr3_keep = self.keep_enriched_clones(CDR_DIC=CDR3_DIC, depth=depth_across_peptide, per_ind_depth=depth_per_individual)
			
			print('NORMALIZING THE FAM NAMES TO BE SAME ACROSS CDR3s')

			norm_alpha_fam_dic = self.norm_fam_cdrs(fam_dic = fam_dic_alpha)
			norm_beta_fam_dic = self.norm_fam_cdrs(fam_dic = fam_dic_beta)

			print('WRITING TO FILE ')

			self.write_outfile(cdr3dic=CDR3b_DIC, cdr3keep=cdr3b_keep, normfamdicalpha=norm_alpha_fam_dic, normfamdicbeta=norm_beta_fam_dic, chain='CDR3B', CDR3_totals=CDR3_totals)
			self.write_outfile(cdr3dic=CDR3a_DIC, cdr3keep=cdr3a_keep, normfamdicalpha=norm_alpha_fam_dic, normfamdicbeta=norm_beta_fam_dic, chain='CDR3A', CDR3_totals=CDR3_totals)
			self.write_outfile(cdr3dic=CDR3_DIC, cdr3keep=cdr3_keep, normfamdicalpha=norm_alpha_fam_dic, normfamdicbeta=norm_beta_fam_dic, chain='CDR3AB', CDR3_totals=CDR3_totals)


if __name__ == '__main__':
	seq_object = sseq_calls(name='SSEQ', filein = filein)
	seq_object.conclude_process(depth_across_peptide=3, depth_per_individual=1)

