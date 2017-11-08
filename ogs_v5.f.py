import re
import os
import subprocess
from sys import argv
import logging
from gff3 import Gff3
import sys

__version__ = '5.0.0'
'''
usage: python ogs_v5.py biosample_ID file.gff3
example: python ogs_v5.py GCA_000696205.1 clec_OGS_v1_2_with_pep_CDS.gff3
'''
#2. remodel pseudogene
def Remodel_pseudogenes(gff_data, logging):

    pseudo_list=[]

    for line in gff_data.lines:
        if line["line_type"] == "feature":
            if "pseudogene" in line["type"]:
                line['type'] = "gene"
                # get gene_ID
                try:
                    gene_ID = line['attributes']['ID']
                    pseudo_list.append(gene_ID)
                except:
                    logging.warning('[Missing ID] - Line: %s' % line['line_index'])


    #for each child of pseudogene, alter col9
    for pi in range(len(pseudo_list)):
        for line in gff_data.lines:
            if line['line_type'] == 'feature':
                try:
                    for attribute_key, attribute_value in line['attributes'].items():
                        if pseudo_list[pi] in str(attribute_value):
                            line['attributes']['pseudogene'] = 'unknown'
                            ## check it's child col3 whether has multi-name , then throw out warning !
                            if (line['type'] != 'gene') and not ( re.search('pseudogenic_transcript', line['type']) or re.search('pseudogenic_exon', line['type'])):
                                #get child_ID
                                try:
                                    child_ID = line['attributes']['ID']
                                    feature = line['type']
                                    logging.warning(" multiple types of child features : "+"ID= "+child_ID+";  feature= "+feature)
                                except:
                                    logging.warning('[Missing ID] - Line: %s' % line['line_index'])
                except:
                    logging.warning('[Missing Attributes] - Line: %s' % line['line_index'])

    logging.info('remodel pseudogenes done!')

#3-1 util check for locus tag in col9
def check_locus(gff_data, f_out, logging):
    idx0 = []
    idx1 = []
    for line in gff_data.lines:
        if line["line_type"] == "feature":
            if "gene" in line["type"]:
                idx0.append(line)
                try:
                    for attribute_key, attribute_value in line['attributes'].items():
                        if "locus_tag" in  attribute_key:
                            idx1.append(line)
                except:
                    logging.warning('[Missing Attributes] - Line: %s' % line['line_index'])

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
def locus_tag_util(biosample, logging):

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
def ID_attribute(locus_tag, gff_data):
    idx1 = []
    for line in gff_data.lines:
        if line['line_type'] == 'feature':
            if "mRNA" in line['type']:
                try:
                    ID = line['attributes']['ID']
                    line['attributes']['transcript_id'] = "gnl|"+locus_tag+"|"+ID
                except:
                    print "mRNA_ID ERROR"
                    print line['line_raw']


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
def Scan_for_genome_alterations(gff_data, f_in, f_out, logging):
    dict={}
    alteration_list=["stop_codon_read_through", "stop_codon_readthrough", "deletion", "insertion", "substitution"]

    for i in range(len(alteration_list)):
        dict[alteration_list[i]] = 0
        for line in gff_data.lines:
            if line['line_type'] == 'feature':
                if alteration_list[i] in line['type']:
                    dict[alteration_list[i]] += 1

  #  print "Program Summary:"
    f_out.write("\nProgram Summary:")
    sum= " ".join([str(v)+" of type "+k+","  for k,v in dict.items()])
   # print "Found the following features in file [input gff3 file]: "+sum.rstrip(",")+"."
    f_out.write( "\nFound the following features in file "+f_in+": "+sum.rstrip(",")+".")
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
def Dbxref(gff_data, logging):
    for line in gff_data.lines:
        if line['line_type'] == 'feature':
            if "gene" in line['type']:
                try:
                    gene_ID = line['attributes']['ID']
                    print gene_ID
                    if not line['attributes'].has_key('Dbxref'):
                        line['attribues']['Dbxref'] = "I5KNAL:" + gene_ID
                except:
                    print "error"
                    logging.warning("Dbxref error gene feature : "+line['line_raw'])
            else:
                # second, for non gene feature, remove Dbxref tag
                try:
                    if line['attributes'].has_key('Dbxref'):
                        del line['attributes']['Dbxref']
                except:
                    logging.warning('[Missing Attribute] - Line %s' % line['line_index'])



def Remodel_pseudogenes_mod(gff_data, logging):

    feature_type = ['pseudogene', 'pseudogenic_transcript','pseudogenic_exon']  # total: 14214, only 2 pseudogenic_transcript
    for line in gff_data.lines:
        if line['line_type'] == 'feature':
            if any(x in line['type'] for x in feature_type):
                try:
                    line['attributes']['pseudogene'] = 'unknown'
                    if (line['type'] != 'pseudogene') and not (re.search('pseudogenic_transcript', line['type']) or re.search('pseudogenic_exon', line['type'])):
                        try:
                            child_ID = line['attributes']['ID']
                            feature = line['type']
                            logging.warning( " multiple types of child features : " + "ID= " + child_ID + "; features= " + feature)
                        except:
                            logging.warning('[Missing Attribute] - Line %s' % line['line_index'])
                    if (line['type'] == 'pseudogene'):
                        line['type'] = 'gene'
                except:
                    logging.warning('[Missing Attribute] - Line %s' % line['line_index'])
    logging.info('remodel pseudogenes done!')



def cds_id_attribute_mod(locus_tag, gff_data, f_out, logging):
        # TODO: if equal to previous ID then pass  checking
        logging.info("start cds_id_attribute")
        df2 = []
        df_pep = []
        df_pepParent_dict = {}
        for line in gff_data.lines:
            if line['line_type'] == 'feature':
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
                        i['attributes']['protein_id'] = "gnl|" + locus_tag + "|" + poly_ID
                        #elif(len( df_pep.where(df_pep[8].str.contains("Parent=" + ID)).dropna().index.tolist())>0):
                elif ID in df_pepParent_dict:
                    if len(df_pepParent_dict[ID]) > 0:
                        print "len a child of cds : "+str(len(df_pepParent_dict[ID]))
                        #poly_df = df_pep.where(df_pep[8].str.contains("Parent=" + ID))
                        #poly_idx = poly_df.dropna().index.tolist()
                        for j in df_pepParent_dict[ID]:
                            poly_ID = j['attributes']['ID']
                            print "poly ID "+ poly_ID
                        i['attributes']['protein_id'] = "gnl|" + locus_tag + "|" + poly_ID
                else:
                    print " polypeptide is niether sharing parent with CDS nor a child of CDS!"
                    logging.warning("polypeptide is niether sharing parent with CDS nor a child of CDS!" )

                    #polyID="".join(re.match("^.*ID=([^;]+).*$", gff_data.iloc[i, 8]).groups())

            except:
                print "CDS_ID ERROR"
                f_out.write("\nCDS_ID error\n")
                f_out.write(i['line_raw'])
                logging.warning("CDS  id error" + i['line_raw'])#2

if __name__ == '__main__':
    import sys
    import argparse
    from textwrap import dedent
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\
    clean Official Gene Sets prior to NCBI submission

    Example:
        python  ogs_v5.py  assembly_ID  gff3_file
    """))

    parser.add_argument('assembly_ID', type=str, help='assembly ID')
    parser.add_argument('gff3_file', type=str, help='input gff3 file')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    args = parser.parse_args()
    base = os.path.basename(args.gff3_file)
    file_name1=os.path.splitext(base)[0]
    file_name2=os.path.splitext(base)[1]
    print args.assembly_ID
    print args.gff3_file

    logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                                            datefmt='%m/%d/%Y %I:%M:%S =%p',
                                            filename=file_name1+'.log',
                                            filemode='w',
                                            level=logging.INFO)
    logging.info('start program')
    f_out=open(file_name1+'_out.txt', 'w')
    gff_data = Gff3(gff_file=args.gff3_file)

    check_locus(gff_data=gff_data, f_out=f_out, logging=logging)#3
    LOCUSTAG=locus_tag_util(biosample=args.assembly_ID, logging=logging)#3
    logging.info('Finish locus tag')
    Remodel_pseudogenes_mod(gff_data=gff_data, logging=logging)#2
    logging.info('Finish remodel pseudogenes')

    Scan_for_genome_alterations(gff_data=gff_data, f_in=args.gff3_file, f_out=f_out, logging=logging) #9
    logging.info('Finish scan for genome alterations')
    ID_attribute(LOCUSTAG, gff_data=gff_data) #4
    logging.info('Finish Id attribute')
    Dbxref(gff_data=gff_data, logging=logging) #12
    logging.info('Finish Dbxref attribute')

    cds_id_attribute_mod(locus_tag=LOCUSTAG, gff_data=gff_data, f_out=f_out, logging=logging) #4
    logging.info('Finish cds_id attribute')
    logging.info('Finish program, start generating new gff3 file.')


    file_out=file_name1+".out"+file_name2

    gff_data.write(file_out)
