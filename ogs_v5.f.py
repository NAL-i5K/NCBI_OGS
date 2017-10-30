import pandas as pd
import re
import os
import subprocess
from sys import argv
import logging
from gff3 import Gff3
import sys
'''
usage: python ogs_v5.py biosample_ID file.gff3
example: python ogs_v5.py GCA_000696205.1 clec_OGS_v1_2_with_pep_CDS.gff3
'''
Biosample=argv[1]
file_in=argv[2]

base=os.path.basename(file_in)
#print base
file_name1=os.path.splitext(base)[0]
file_name2=os.path.splitext(base)[1]

#file_in='clec_OGS_v1_2_with_pep_CDS.gff3'
print Biosample
print file_in

logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    datefmt='%m/%d/%Y %I:%M:%S =%p',
                    filename=file_name1+'.log',
                    filemode='w',
                    level=logging.INFO)
logging.info('start program')
#gff_data = pd.read_csv('clec_OGS_v1_2_with_pep_CDS.gff3', sep="\t", header = None,comment='#')
'''
remove NOTES attribute
name attributes: gene, mRNA, pseudogenic_transcript

'''

f_out=open(file_name1+'_out.txt', 'w')
#1. remove %09
#preprocess= subprocess.Popen("sed s/%09//g "+ file_in +">temp.gff", stdout=subprocess.PIPE,shell=True)
#gff_data = pd.read_csv('temp.gff', sep="\t", header = None,comment='#')
gff_data = pd.read_csv(file_in, sep="\t", header = None,comment='#')
gff = Gff3(gff_file=file_in)
#2. remodel pseudogene
def Remodel_pseudogenes():
    pseudo_list=[]


    df0 = gff_data.where(gff_data[2].str.contains("pseudogene"))
    idx0 = df0.dropna().index.tolist()
    for index in idx0:
        gff_data.iloc[index, 2] = "gene"
       # gff_data.iloc[index, 8] = gff_data.loc[index, 8] + ";pseudogene=unknown"
        # get gene_ID
        gene_ID = "".join(re.match("^.*ID=([^;]+);.+$", gff_data.iloc[index, 8]).groups())
        # print gene_ID
        pseudo_list.append(gene_ID)

    #pseudo_list=['CLEC026002', 'CLEC026001','CLEC0260214']

    #for each child of pseudogene, alter col9
    for pi in range(len(pseudo_list)):
        df1 = gff_data.where(gff_data[8].str.contains(pseudo_list[pi]))
        idx2 = df1.dropna().index.tolist()
        for i in idx2:
            gff_data.iloc[i, 8] = str(gff_data.loc[i, 8]).rstrip(';') + ";pseudogene=unknown"
        ## check it's child col3 whether has multi-name , then throw out warning !
           # if (gff_data.loc[i,2]!="gene") and (gff_data.loc[i,2]!="pseudogenic_transcript" or gff_data.loc[i,2]!= "pseudogenic_exon"):
            if (gff_data.loc[i, 2] != "gene") and not ( re.search( "pseudogenic_transcript",gff_data.loc[i, 2]) or re.search( "pseudogenic_exon",gff_data.loc[i, 2])):
                child_ID = "".join(re.match("^.*ID=([^;]+);.+$", gff_data.iloc[i, 8]).groups())
                feature=gff_data.loc[i,2]
                logging.warning(" multiple types of child features : "+"ID= "+child_ID+";  feature= "+feature)

    logging.info('remodel pseudogenes done!')
   # gff_data.to_csv('clec_OGS_v1_2.step2.gff3', header=False, index=False, sep='\t')
#3-1 util check for locus tag in col9
def check_locus():
    df0 = gff_data.where(gff_data[2].str.contains("gene"))
    idx0 = df0.dropna().index.tolist()
    df1 = df0.where(df0[8].str.contains("locus_tag"))
    idx1 = df1.dropna().index.tolist()
    if len(idx1) == len(idx0):
        f_out.write("All gene features contain a locus_tag attribute")
        logging.info("All gene features contain a locus_tag attribute")
    elif 0 < len(idx1) < len(idx0):
        f_out.write("Some gene features contain a locus_tag attribute")
        logging.info("Some gene features contain a locus_tag attribute")
    elif len(idx1) == 0:
        f_out.write("No gene features contain a locus_tag attribute")
        logging.info("No gene features contain a locus_tag attribute")
    logging.info('check locus done!')
#3-2. util get locus tag
def locus_tag_util(biosample):

    import subprocess
    import logging
    import re

    biosample_1 = subprocess.Popen("esearch -db assembly -query "+biosample, stdout=subprocess.PIPE, shell=True)
    biosample_2 = subprocess.Popen("efetch -format docsum", stdin=biosample_1.stdout, stdout=subprocess.PIPE,
                                   shell=True)
    # biosample_3=subprocess.Popen("xtract -pattern DocumentSummary -element BiosampleID", stdin=biosample_2.stdout,stdout=subprocess.PIPE, shell=True)
    biosample_3 = subprocess.Popen("xtract -pattern DocumentSummary -element BioSampleAccn", stdin=biosample_2.stdout,
                                   stdout=subprocess.PIPE, shell=True)
    proc_stdout_biosample = biosample_3.communicate()[0].strip()


    process = subprocess.Popen("esearch -db assembly -query " + biosample ,stdout=subprocess.PIPE, shell=True)
    process2= subprocess.Popen("efetch -format docsum" ,stdin=process.stdout , stdout=subprocess.PIPE,shell=True)
    process3=subprocess.Popen("xtract -pattern DocumentSummary -element BioprojectId", stdin=process2.stdout,stdout=subprocess.PIPE, shell=True)
    proc_stdout_bioproject = process3.communicate()[0].strip()
#167479  348318
   # print "bioproject\n"+ proc_stdout_bioproject
    bioprojectid_list=proc_stdout_bioproject.split("\t")
   # print bioprojectid_list
    locus_dict = {}
    for bioproject_i in bioprojectid_list:
        process31=subprocess.Popen("efetch -db bioproject -id "+bioproject_i+" -format xml",stdout=subprocess.PIPE, shell=True)
        #process32=subprocess.Popen("xtract -pattern DocumentSummary -element LocusTagPrefix", stdin=process31.stdout,stdout=subprocess.PIPE, shell=True)

        process32=subprocess.Popen("grep LocusTagPrefix -", stdin=process31.stdout,stdout=subprocess.PIPE, shell=True)

        proc_stdout_locus = process32.communicate()[0].strip()
       # print "locus: "+proc_stdout_locus
        locus_list = proc_stdout_locus.lstrip("\t").split("\n")
       # print "locus_list: "+str(locus_list)
        for i in range(len(locus_list)):
            if(re.match(".*<LocusTagPrefix.*biosample_id=\"(\w+)\".+>", locus_list[i])):
                biosample_id = "".join(re.match(".*<LocusTagPrefix.*biosample_id=\"(\w+)\".+>", locus_list[i]).groups())
                locustag = "".join(re.match(".*<LocusTagPrefix.*biosample_id=\"\w+\">\s*([\w]+).+>", locus_list[i]).groups())
       #         print biosample_id,locustag
                locus_dict.update({biosample_id: locustag})
            elif(re.match(".*<LocusTagPrefix>(\w+)</LocusTagPrefix>", locus_list[i])):
                locustag="".join(re.match(".*<LocusTagPrefix>(\w+)</LocusTagPrefix>", locus_list[i]).groups())
                locus_dict.update({proc_stdout_biosample:locustag})
    try:
        LOCUS_TAG = locus_dict.get(proc_stdout_biosample)
      #  print LOCUS_TAG
        logging.info("Locus tag: "+LOCUS_TAG)
        return LOCUS_TAG

    except:
        logging.warning("biosample ID not in locus prefix list")
#4. ID attribute  mRNAS , CDS transcript_id, protein_id
def ID_attribute(locus_tag):
    df1 = gff_data.where(gff_data[2].str.contains("mRNA"))
    idx1 = df1.dropna().index.tolist()
    for i in idx1:
        try:
            ID="".join(re.match("^.*ID=([^;]+).*$", gff_data.iloc[i, 8]).groups())
            gff_data.iloc[i, 8] = str(gff_data.loc[i, 8]).rstrip(";") + ";transcript_id=gnl|"+locus_tag+"|"+ID
            #print gff_data.iloc[i,8]
        except:
            print "mRNA_ID ERROR"
            print  gff_data.iloc[i,8]



#5. remove note attribute
def remove_Note():
        df1 = gff_data.where(gff_data[8].str.contains("Note="))
        idx2 = df1.dropna().index.tolist()
        ## check if the Note is the last term!
        for i in idx2:
            try:
                Note_term = "".join(re.match("^.*Note=([^;]+);.*$", gff_data.iloc[i, 8]).groups())
              #  print Note_term
                Note_new= re.sub(re.escape(Note_term), "Additional notes are available on the i5k Workspace gene page.", gff_data.iloc[i,8])
              #  print Note_new
                gff_data.iloc[i, 8] = Note_new
              #  print gff_data.iloc[i,8]
            except:
                print gff_data.iloc[i,8]
                print "".join(re.match("^.*Note=([^;]+).*$", gff_data.iloc[i, 8]).groups())
                logging.warning('remove Note error: '+ gff_data.iloc[i,8])
#		print "".join(re.match("^.*Note=([^;]+).*$", gff_data.iloc[i, 8]).groups())

#6. Name attributes(gene)
def name_Attributes_gene():
        df1 = gff_data.where(gff_data[2].str.contains("gene"))
        idx1 = df1.dropna().index.tolist()
        #print len(idx1) # 14087
        df2 = df1.where(df1[8].str.contains("Name="))
        idx2 = df2.dropna().index.tolist()
        #print idx2,len(idx2) #696
        for i in idx2:
                try:
                        name_attribute = "".join(re.match("^.*Name=([^;]+).*$", gff_data.iloc[i, 8]).groups()).rstrip("\s")
                        ## TODO:
    #                    print name_attribute
                        temp="Name="+name_attribute+";gene="+name_attribute
    #                    print temp
                        name_new = re.sub("Name="+name_attribute, (temp), gff_data.iloc[i, 8])
    #                    print name_new
                        gff_data.iloc[i, 8] = name_new
                except:
                        logging.warning( "ERROR gene name attribute:"+gff_data.iloc[i,8])
                        name_attribute= "".join(re.match("^.*Name=([^;]+).*$", gff_data.iloc[i, 8]).groups()).rstrip("\s")
                        new = name_attribute+";gene="+name_attribute
                        #a= re.sub(name_attribute, new, gff_data.iloc[i, 8])
                        print new

#7. Name_attribute for mRNA, pseudo_transcript
# multi conditions for features:
def name_Attribute_mRNA_pseudo():
        feature_type = ['mRNA', 'pseudogenic_transcript']  #total: 14214, only 2 pseudogenic_transcript
        df1 = gff_data.where(gff_data[2].str.contains('|'.join(feature_type)))
        idx1 = df1.dropna().index.tolist()
        df2 = df1.where(df1[8].str.contains("Name="))
        idx2 = df2.dropna().index.tolist()
        #print df2, len(idx2) # 1481
        for i in idx2:
                 try:
                        name_attribute = "".join(re.match("^.*Name=([^;]+).*$", gff_data.iloc[i, 8]).groups())
                        temp="Name="+name_attribute+";product="+name_attribute
                        name_new = re.sub("Name="+name_attribute, (temp), gff_data.iloc[i, 8])
                        #print name_new
                        gff_data.iloc[i, 8] = name_new
                 except:
                        if (gff_data.iloc[i,2]=="mRNA"):
                            logging.warning( "ERROR mRNA name attribute:"+gff_data.iloc[i,8])
                            print "".join(re.match("^.*Name=([^;]+).*$", gff_data.iloc[i, 8]).groups())

                        elif (gff_data.iloc[i,2]=="pseudogenic_transcript"):

                            logging.warning( "ERROR pseudogenic_transcript name attribute:"+gff_data.iloc[i,8])
                            print "".join(re.match("^.*Name=([^;]+).*$", gff_data.iloc[i, 8]).groups())
#8. name attribute (tRNAs, rRNAs and ncRNAs)
def name_Attribute_trncRNA():
    feature_type = ['tRNA', 'rRNA','ncRNA']  #total: 14214, only 2 pseudogenic_transcript
    df1 = gff_data.where(gff_data[2].str.contains('|'.join(feature_type)))
    idx1 = df1.dropna().index.tolist()
    #print df1 #3
    df2 = df1.where(df1[8].str.contains("Name="))
    idx2 = df2.dropna().index.tolist()
    for i in idx2:
        try:
            name_attribute = "".join(re.match("^.*Name=([^;]+).+$", gff_data.iloc[i, 8]).groups())
            name_new = re.sub("Name=", "product=" , gff_data.iloc[i, 8])
           # print name_new
            gff_data.iloc[i, 8] = name_new
        except:

            logging.warning("ERROR rRNA,tRNA,ncRNA  name attribute:" + gff_data.iloc[i, 8])
            print "".join(re.match("^.+Name=([^;]+).+$", gff_data.iloc[i, 8]).groups())

#9. stop_codon_read_through, stop_codon_readthrough, deletion, insertion, substitution
def Scan_for_genome_alterations():
    dict={}
    alteration_list=["stop_codon_read_through", "stop_codon_readthrough", "deletion", "insertion", "substitution"]
    for i in range(len(alteration_list)):
        df1 = gff_data.where(gff_data[2].str.contains(alteration_list[i]))
        dx1 = df1.dropna().index.tolist()
        dict.update({alteration_list[i]:len(dx1)})
  #  print "Program Summary:"
    f_out.write("\nProgram Summary:")
    sum= " ".join([str(v)+" of type "+k+","  for k,v in dict.items()])
   # print "Found the following features in file [input gff3 file]: "+sum.rstrip(",")+"."
    f_out.write( "\nFound the following features in file"+file_in+": "+sum.rstrip(",")+".")
    #print "You may need to manually add a transl_except attribute to the corresponding CDS feature lines. See https://www.ncbi.nlm.nih.gov/sites/genbank/genomes_gff/ for details."
    f_out.write("\nYou may need to manually add a transl_except attribute to the corresponding CDS feature lines. See https://www.ncbi.nlm.nih.gov/sites/genbank/genomes_gff/ for details.")

#10. Description attributes
def description_attribute():
    df1 = gff_data.where(gff_data[8].str.contains("description="))
    dx1 = df1.dropna().index.tolist() #1254
    df2= gff_data.where(gff_data[2].str.contains("gene"))
    df21=df2.where(df2[8].str.contains("description="))
    dx2 = df21.dropna().index.tolist() #425
    #print len(dx1), len(dx2)
    #print str(len(dx1)-len(dx2))+"  description attributes found outside of gene features."
    f_out.write("\n"+str(len(dx1)-len(dx2))+"  description attributes found outside of gene features.")

#11. Symbol
#Program will scan for symbol attributes in features other than gene features.
def symbol_attribute():
    df1 = gff_data.where(gff_data[8].str.contains("symbol="))
    dx1 = df1.dropna().index.tolist() #1254

    df2= gff_data.where(gff_data[2].str.contains("gene"))
    df21=df2.where(df2[8].str.contains("description="))
    dx2 = df21.dropna().index.tolist() #425
    #print len(dx1), len(dx2)
    print str(len(dx1)-len(dx2))+"  symbol attributes found outside of gene features."
    f_out.write("\n"+str(len(dx1)-len(dx2))+"  symbol attributes found outside of gene features.")

#12. Dbxrefs
## first for gene feature, copy ID to Dbxref
def Dbxref():
    df1=gff_data.where(gff_data[2].str.contains("gene"))
    dx1 = df1.dropna().index.tolist() #1254
    for i in dx1:
        try:
            gene_ID = "".join(re.match("^.*ID=([^;]+).*$", gff_data.iloc[i, 8]).groups())
            print gene_ID

            # print name_new
          #  if not (re.search("Dbxref=I5KNAL:",gff_data.iloc[i,8])):
            if (re.match("Dbxref=I5KNAL:", gff_data.iloc[i, 8]) is None ): 
                gff_data.iloc[i, 8] =  gff_data.iloc[i, 8] +";Dbxref:I5KNAL:"+gene_ID

        except:
            print "error"
            logging.warning("Dbxref error gene feature : "+gff_data.iloc[i,8])
    ## second, for non gene feature, remove Dbxref tag
    df2=gff_data.where(gff_data[2].str.contains("gene"))
    dx2 = df2.dropna().index.tolist() #1254
    total=list(range(len(gff_data)))
    ret_list = list(set(dx2)^set(total)) #non gene feature
    for i in ret_list:
        
        if re.search("Dbxref",gff_data.iloc[i,8]):
            try:
                if (re.match("^.*(Dbxref=[^;]+);.*$", gff_data.iloc[i, 8]) ):
                    Dbxref= "".join(re.match("^.*(Dbxref=[^;]+).*$", gff_data.iloc[i, 8]).groups())+";"
                    new_dbxref=re.sub(re.escape(Dbxref),"",gff_data.iloc[i,8])
                    gff_data.iloc[i,8]=new_dbxref
                   # print gff_data.iloc[i,8]
                   # print new_dbxref
                elif (re.match("^.*(Dbxref=[^;]);$", gff_data.iloc[i, 8])):
                    Dbxref = "".join(re.match("^.*(Dbxref=[^;]+);$", gff_data.iloc[i, 8]).groups()) + ";"
                    new_dbxref = re.sub(re.escape(Dbxref), "", gff_data.iloc[i, 8])
                    gff_data.iloc[i, 8] = new_dbxref
                elif(re.match("^.*(Dbxref=[^;])$", gff_data.iloc[i, 8])):
                    Dbxref = "".join(re.match("^.*(Dbxref=[^;]+)$", gff_data.iloc[i, 8]).groups())
                    new_dbxref = re.sub(re.escape(Dbxref), "", gff_data.iloc[i, 8])
                    gff_data.iloc[i, 8] = new_dbxref
            except:
                    logging.warning("Dbxref error non gene feature remove Dbxref tag: "+ gff_data.iloc[i,8] )
        

def Remodel_pseudogenes_mod():

    feature_type = ['pseudogene', 'pseudogenic_transcript','pseudogenic_exon']  # total: 14214, only 2 pseudogenic_transcript
    df_mo = gff_data.where(gff_data[2].str.contains('|'.join(feature_type)))
    idx_mo = df_mo.dropna().index.tolist()
    for i in idx_mo:
        gff_data.iloc[i, 8] = str(gff_data.loc[i, 8]).rstrip(';') + ";pseudogene=unknown"
        if (gff_data.loc[i, 2] != "pseudogene") and not (re.search("pseudogenic_transcript", gff_data.loc[i, 2]) or re.search("pseudogenic_exon", gff_data.loc[i, 2])):
            child_ID = "".join(re.match("^.*ID=([^;]+);.+$", gff_data.iloc[i, 8]).groups())
            feature = gff_data.loc[i, 2]
            logging.warning(" multiple types of child features : " + "ID= " + child_ID + ";  feature= " + feature)
        if (gff_data.loc[i, 2] == "pseudogene"):
            gff_data.iloc[i, 2] = "gene"

    logging.info('remodel pseudogenes done!')


def cds_id_attribute_mod(locus_tag):
        # TODO: if equal to previous ID then pass  checking
        logging.info("start cds_id_attribute")
        df2 = []
        df_pep = []
        df_pepParent_dict = {}
        for line in gff.lines:
            if "CDS" in line['type']:
                df2.append(line)
            if "polypeptide" in line['type']:
                df_pep.append(line)
                try:
                    if ",".join(line['attributes']['Parent']) not in df_pepParent_dict:
                        df_pepParent_dict[",".join(line['attributes']['Parent'])] = [line]
                    else:
                        df_pepParent_dict[",".join(line['attributes']['Parent'])].append(line)
                except:
                    pass
  
        #df2 = gff_data.where(gff_data[2].str.contains("CDS"))
        #df2 = [line for line in gff.lines if "CDS" in line['type']]
        #idx2 = df2.dropna().index.tolist()
        #df_pep = gff_data.where(gff_data[2].str.contains("polypeptide"))
        #df_pep = [line for line in gff.lines if "polypeptide" in line['type']] 
        for i in df2:
                try:
                    #Parent_ID_form = "".join(re.match("^.*Parent=([^;]+).*$", gff_data.iloc[i, 8]).groups())
                    Parent_ID_form = ",".join(i['attributes']['Parent'])
                    print "Parent : "+  Parent_ID_form
                    #ID = "".join(re.match("^.*ID=([^;]+).*$", gff_data.iloc[i, 8]).groups())  # Parent=CLEC000001-RA
                    ID = i['attributes']['ID']
                    print "ID: "+ID
                    poly_ID =""
                    
                    #if (len(df_pep.where(df_pep[8].str.contains("Parent=" + Parent_ID_form)).dropna().index.tolist()) > 0):
                    if Parent_ID_form in df_pepParent_dict:
                        if len(df_pepParent_dict[Parent_ID_form]) > 0:
                            print "len share p: "+str(len(df_pepParent_dict[Parent_ID_form]))
                        
                            #poly_df=df_pep.where(df_pep[8].str.contains("Parent=" + Parent_ID_form))
                            #poly_idx = poly_df.dropna().index.tolist()
                            for j in df_pepParent_dict[Parent_ID_form]:
                               poly_ID= j['attributes']['ID']
                               print "poly ID "+ poly_ID
                            gff_data.iloc[i['line_index']-1, 8] = str(gff_data.loc[i['line_index']-1, 8]).rstrip(";") + ";protein_id=gnl|" + locus_tag + "|" + poly_ID
                        #elif(len( df_pep.where(df_pep[8].str.contains("Parent=" + ID)).dropna().index.tolist())>0):
                    elif ID in df_pepParent_dict:
                        if len(df_pepParent_dict[ID]) > 0:
                            print "len a child of cds : "+str(len(df_pepParent_dict[ID]))
                            #poly_df = df_pep.where(df_pep[8].str.contains("Parent=" + ID))
                            #poly_idx = poly_df.dropna().index.tolist()
                            for j in df_pepParent_dict[ID]:
                                poly_ID = j['attributes']['ID']
                                print "poly ID "+ poly_ID
                            gff_data.iloc[i['line_index']-1, 8] = str(gff_data.loc[i['line_index']-1, 8]).rstrip(";") + ";protein_id=gnl|" + locus_tag + "|" + poly_ID
                    else:
                        print " polypeptide is niether sharing parent with CDS nor a child of CDS!"
                        logging.warning("polypeptide is niether sharing parent with CDS nor a child of CDS!" )
    
                           #polyID="".join(re.match("^.*ID=([^;]+).*$", gff_data.iloc[i, 8]).groups())

                except:
                    print "CDS_ID ERROR"
                    f_out.write("\nCDS_ID error\n")
                    f_out.write(gff_data.iloc[i['line_index']-1, 8])
                    logging.warning("CDS  id error" + gff_data.iloc[i['line_index']-1,8])#2
check_locus()#3
LOCUSTAG=locus_tag_util(Biosample)#3 
logging.info('Finish locus tag')
Remodel_pseudogenes_mod()#2
logging.info('Finish remodel pseudogenes')

#remove_Note() #5
#logging.info('Finish remove Note')
#name_Attributes_gene()#6
#name_Attribute_mRNA_pseudo() #7
#name_Attribute_trncRNA() #8
#logging.info('Finish name attribute')
Scan_for_genome_alterations() #9
logging.info('Finish scan for genome alterations')
#description_attribute() #10
#logging.info('Finish description attribute')
#symbol_attribute() #11
#logging.info('Finish symbol attribute')
ID_attribute(LOCUSTAG) #4
logging.info('Finish Id attribute')
#cds_id_attribute(LOCUSTAG) #4
#logging.info('Finish cds_id attribute')
Dbxref() #12
logging.info('Finish Dbxref attribute')

cds_id_attribute_mod(LOCUSTAG) #4
logging.info('Finish cds_id attribute')
logging.info('Finish program, start generating new gff3 file.')


file_out=file_name1+".out"+file_name2

gff_data.to_csv(file_out, header=False, index=False, sep='\t')


