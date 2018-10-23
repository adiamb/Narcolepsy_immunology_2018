import re
import copy
import gzip
import os
import sys 
import datetime 
import numpy as np
import scipy
#file = '/media/labcomp/HDD2/RCL_FILES_AUG23_2017/SINGLE_CELL_SEQ_COMP_GUO/FAMILY_USAGE_SSEQ/RAW_HCRT1_8_HCRT2_6_TO_FIGURE.csv'
file = '../COMBINED_SSEQ_GUO_October12_2018_edit_fornorm_lowthreshold.csv'
dttime=datetime.datetime.now().strftime ("%Y%m%d")
print dttime

#### define subroutines
def make_CDRDIC(indic, makekey, cdr):
	if cdr in indic:
		get_exist_cdr = indic.get(cdr)
		if makekey in get_exist_cdr:
			get_count = get_exist_cdr.get(makekey)
			indic[cdr][makekey] = get_count +1
		else:
			indic[cdr][makekey] = 1
	else:
		indic[cdr]={makekey:1}

## 
def counts_dic(makekey, total_dic, comb_call):
	"""this function counts the totals of valid CDR3s both a and b and combined returns a dic that has totals"""
	if makekey in total_dic:
		get_call = total_dic.get(make_key)
		if comb_call in get_call:
			get_count = get_call.get(comb_call)
			total_dic[make_key][comb_call] = get_count + 1
		else:
			total_dic[make_key][comb_call] =  1
	else:
		total_dic[make_key] = {comb_call:1}


### 
def family_specific(make_key, indic, fam_to_be_processed, cdr):
	"""ths fucntion compliles  family specific analysis"""
	if make_key in indic:
		get_exist_cdr = indic.get(make_key)
		if cdr in get_exist_cdr:
			get_cdr_count = get_exist_cdr.get(cdr)
			indic[make_key][cdr] = get_cdr_count+1
		else:
			indic[make_key][cdr] = 1
	else:
		indic[make_key]= {cdr:1}

##
def norm_fam_cdrs(fam_dic):
	"""will take a cdr3 and associated families and normalize the fam names for same cdrs based on the counts"""
	norm_cdrs ={}
	for k, v in fam_dic.items():
		for k2, v2 in v.items():
			if len(v) > 1:
				print k, k2, v2
				if v2 > 1:
					norm_cdrs[k] = k2
				else: ##### incase same count assig the first instance of the fam to main call
					if k not in norm_cdrs:
						norm_cdrs[k] = k2
			else:
				norm_cdrs[k] = k2
	return norm_cdrs

##
def keep_enriched_clones(CDR_DIC, depth, per_ind_depth):
	"""# takes a master dictionary and processses with a condition that a cdr3 must be present atleast in 
	3 individuals and more than 3x returns a keep list of cdrs"""
	out_cdr={}
	for k, v in CDR_DIC.items():
		temp_dic ={}
		if len(v) > per_ind_depth:
			for k2, v2 in v.items():
				parse_k2 =  k2.split(',')
				dbid = parse_k2[0]+','+parse_k2[-1]
				pep = parse_k2[1]
				if k in out_cdr:
					get_cdr_dic = out_cdr.get(k)
					if pep in get_cdr_dic:
						get_pep_dic = get_cdr_dic.get(pep)
						if dbid in get_pep_dic:
							get_count = get_pep_dic.get(dbid)
							out_cdr[k][pep][dbid] = get_count+v2
						else:
							out_cdr[k][pep][dbid] = v2
					else:
						out_cdr[k][pep] = {dbid:v2}
				else:
					out_cdr[k]= {pep:{dbid:v2}}
	keep_list = []					
	for k, v in out_cdr.items():
		for k2, v2 in v.items():
			if len(v2) >= depth:
				print k, k2
				if k not in keep_list:
					keep_list.append(k)
	return keep_list

## 
def get_uniques_per_cdr3(CDR, CDR_DIC):
	"""this function writes out all the unique CDR3s found for each peptide"""
	uniques_DIC={}
	dttime=datetime.datetime.now().strftime ("%Y%m%d")
	for k, v in CDR_DIC.items():
		for k2 in v.iterkeys():
			make_key=k2.split(',')[1]
			if 'Blank' not in make_key:
				if 'NP23' not in make_key:
					if make_key in uniques_DIC:
						get_count = uniques_DIC.get(make_key)
						uniques_DIC[make_key] = get_count +1
					else:
						uniques_DIC[make_key] =1

	makefile = 'CDR3_SSEQ_UNIQUES_TOTAL_'+dttime+'.csv'
	if os.path.exists(makefile):
		outfile = open(makefile, 'a')
		
	else:
		outfile = open(makefile, 'w')
		outfile.write('PEPTIDE,COUNTS_UNIQUE,CDR'+'\n')
	for k, v in uniques_DIC.items():
		outfile.write(k+','+str(v)+','+CDR+'\n')
	outfile.close()

##
def get_normalized_fam(cdr, normfamdicalpha, normfamdicbeta, chain):
	"""fucntion returns a normalized familiy of a given cdr3"""
	#global make_k
	if chain == 'CDR3AB':
		get_alpha_CDR = cdr.split(':')[0]
		get_beta_CDR = cdr.split(':')[1]
		if get_alpha_CDR in normfamdicalpha:
			get_fam = normfamdicalpha.get(get_alpha_CDR).split('_')
			make_k1 = get_fam[0]+'_'+get_alpha_CDR+'_'+get_fam[1]
		if get_beta_CDR in normfamdicbeta:
			get_fam = normfamdicbeta.get(get_beta_CDR)
			make_k2 = get_fam+'_'+get_beta_CDR
		make_k = make_k1+':'+make_k2
	elif chain == 'CDR3A':
		if cdr in normfamdicalpha:
			get_fam = normfamdicalpha.get(cdr).split('_')
			make_k = get_fam[0]+'_'+cdr+'_'+get_fam[1]
	else:
		if cdr in normfamdicbeta:
			get_fam = normfamdicbeta.get(cdr)
			make_k = get_fam+'_'+cdr
	return make_k



##
def write_outfile(cdr3dic, cdr3keep, normfamdicalpha, normfamdicbeta, chain, CDR3_totals):
	""" this function writes out the sharing results """
	dttime=datetime.datetime.now().strftime ("%Y%m%d")
	header=',DBID,PEPTIDE,DX,COUNTS,Totalcalls'+'\n'
	if chain == 'CDR3AB':
		outfile = open('CDR3AB_SHARING_SSEQ_GUO_'+dttime+'.csv', 'w')
		outfile.write('CDR3AB')
		cdr3_total_key = 'both'
	elif chain =='CDR3A':
		outfile = open('CDR3A_SHARING_SSEQ_GUO_'+dttime+'.csv', 'w')
		outfile.write('CDR3A')
		cdr3_total_key = 'a'
	else:
		outfile = open('CDR3B_SHARING_SSEQ_GUO_'+dttime+'.csv', 'w')
		outfile.write('CDR3B')
		cdr3_total_key = 'b'
	outfile.write(header)
	for k, v in cdr3dic.items():
		if k in cdr3keep:
			make_k=get_normalized_fam(cdr = k, normfamdicbeta=normfamdicbeta, normfamdicalpha=normfamdicalpha, chain=chain)
			if len(v) > 1:
				for k2, v2 in v.items():
					if k2 in CDR3_totals:
						if cdr3_total_key in CDR3_totals.get(k2):
							get_counts = CDR3_totals.get(k2).get(cdr3_total_key)
						if 'Blank,Blank,Blank' not in k2:
							outfile.write(make_k+','+k2+','+str(v2)+','+str(get_counts)+','+'\n')
	outfile.close()




################################################################### main calls ###########################################################
usage_headers=[]
DBID_list =[]
fam_dic_alpha={}
fam_dic_beta={}
avoid_ =['DbID', 'Blank']
peptide_data ={}
header =''
CDR3_DIC={}
CDR3a_DIC ={}
CDR3b_DIC ={}
CDR3_totals ={}
with open(file) as infile:
	for line in infile:
		if ";" in line:
			parse_line=line.strip('\n').split(';')
		else:
			parse_line=line.strip('\n').split(',')
		#if parse_line[9] or parse_line[15]:
		if parse_line[1] and parse_line[2] and parse_line[3] and parse_line[0]:
			make_key = ','.join([parse_line[1], parse_line[2].strip(), parse_line[3]])
			#make_key = re.sub(r'(HCRT1-8|HCRT2-6)', 'HCRT', make_key)
			print make_key

		
		if 'DbID' not in make_key:
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
				if cdr:
					make_CDRDIC(indic =CDR3_DIC, makekey=make_key, cdr=cdr)
				comb_call = 'both'
				counts_dic(makekey=make_key, comb_call=comb_call, total_dic=CDR3_totals)
			elif cdr3a_alt and cdr3b and not cdr3a:
				# cdr = VAalt+'_'+cdr3a_alt+'_'+JAalt+':'+VB+'_'+cdr3b
				cdr = cdr3a_alt+':'+cdr3b
				comb_call = 'both'
				make_CDRDIC(indic =CDR3_DIC, makekey=make_key, cdr=cdr)
				counts_dic(makekey=make_key, comb_call=comb_call, total_dic=CDR3_totals)
			elif cdr3a_alt and cdr3a and cdr3b:
				# cdr1 = VAalt+'_'+cdr3a_alt+'_'+JAalt+':'+VB+'_'+cdr3b
				# cdr2 = VA+'_'+cdr3a+'_'+JA+':'+VB+'_'+cdr3b
				cdr1 = cdr3a_alt+':'+cdr3b
				cdr2 =cdr3a+':'+cdr3b
				print "FOUND alt alphas ALPHA {}, alt ALPHA {} BETAS {}".format(cdr3a, cdr3a_alt, cdr3b)
				for cdr in [cdr1, cdr2]:
					comb_call = 'both'
					make_CDRDIC(indic =CDR3_DIC, makekey=make_key, cdr=cdr)
					counts_dic(makekey=make_key, comb_call=comb_call, total_dic=CDR3_totals)
			
			if cdr3a and not cdr3a_alt:######################################################### here this script picks up only alpha chains #############
				#cdr3a = cdr3a
				make_CDRDIC(indic =CDR3a_DIC, makekey=make_key, cdr=cdr3a)
				comb_call = 'a'
				counts_dic(makekey=make_key, comb_call=comb_call, total_dic=CDR3_totals)
				make_CDRDIC(indic=fam_dic_alpha, makekey=VA+'_'+JA, cdr=cdr3a)


			elif cdr3a_alt and not cdr3a:
				#cdr3a_alt= VAalt+'_'+cdr3a_alt+'_'+JAalt
				make_CDRDIC(indic =CDR3a_DIC, makekey=make_key, cdr=cdr3a_alt)
				comb_call = 'a'
				counts_dic(makekey=make_key, comb_call=comb_call, total_dic=CDR3_totals)
				make_CDRDIC(indic=fam_dic_alpha, makekey=VAalt+'_'+JAalt, cdr=cdr3a_alt)
				

			elif cdr3a_alt and cdr3a:
				for cdr, V, J in zip([cdr3a_alt, cdr3a], [VAalt, VA], [JAalt, JA]):
					#cdr = V + '_' +cdr+'_'+J
					make_CDRDIC(indic =CDR3a_DIC, makekey=make_key, cdr=cdr)
					comb_call = 'a'
					counts_dic(makekey=make_key, comb_call=comb_call, total_dic=CDR3_totals)
					make_CDRDIC(indic=fam_dic_alpha, makekey=V+'_'+J, cdr=cdr)

			if cdr3b:######################################################### here this script picks up only beta chains #############
				#cdr3b = VB+'_'+cdr3b
				make_CDRDIC(indic =CDR3b_DIC, makekey=make_key, cdr=cdr3b)
				comb_call = 'b'
				counts_dic(makekey=make_key, comb_call=comb_call, total_dic=CDR3_totals)
				make_CDRDIC(indic=fam_dic_beta, makekey=VB, cdr=cdr3b)

############## normlaize the families for each cdr3 found in the database#######################################
norm_alpha_fam_dic = norm_fam_cdrs(fam_dic = fam_dic_alpha)
norm_beta_fam_dic = norm_fam_cdrs(fam_dic = fam_dic_beta) 

outfile = open('NORMALIZED_ALPHA_CDRS_FAMILIES_'+dttime+'.csv', 'w')
for k, v in norm_alpha_fam_dic.items():
	outfile.write(k+','+v.replace('_', ',')+'\n')
outfile.close()

outfile = open('NORMALIZED_BETA_CDRS_FAMILIES_'+dttime+'.csv', 'w')
for k, v in norm_beta_fam_dic.items():
	outfile.write(k+','+v.replace('_', ',')+'\n')
outfile.close()


############ keep only the enriched clones ######################################## the PNAS paper used a depth of 3 i.e a clone must be present in atleast 3 individuals
cdr3a_keep=keep_enriched_clones(CDR_DIC = CDR3a_DIC, depth=3, per_ind_depth=1)
cdr3b_keep=keep_enriched_clones(CDR_DIC = CDR3b_DIC, depth=3, per_ind_depth=1)
cdr3_keep = keep_enriched_clones(CDR_DIC=CDR3_DIC, depth=3, per_ind_depth=1)

############ write out the totals for each peptide
get_uniques_per_cdr3(CDR='CDR3A', CDR_DIC=CDR3a_DIC)
get_uniques_per_cdr3(CDR='CDR3AB', CDR_DIC=CDR3_DIC)
get_uniques_per_cdr3(CDR='CDR3B', CDR_DIC=CDR3b_DIC)

########## write out the results ###########
write_outfile(cdr3dic=CDR3b_DIC, cdr3keep=cdr3b_keep, normfamdicalpha=norm_alpha_fam_dic, normfamdicbeta=norm_beta_fam_dic, chain='CDR3B', CDR3_totals=CDR3_totals)
write_outfile(cdr3dic=CDR3a_DIC, cdr3keep=cdr3a_keep, normfamdicalpha=norm_alpha_fam_dic, normfamdicbeta=norm_beta_fam_dic, chain='CDR3A', CDR3_totals=CDR3_totals)
write_outfile(cdr3dic=CDR3_DIC, cdr3keep=cdr3_keep, normfamdicalpha=norm_alpha_fam_dic, normfamdicbeta=norm_beta_fam_dic, chain='CDR3AB', CDR3_totals=CDR3_totals)
