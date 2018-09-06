## NCBI_OGS
* [gitlab issue](https://gitlab.com/i5k_Workspace/workspace_roadmap/issues/406)
* clean Official Gene Sets prior to NCBI submission  
* Main contributors are [Chao-I Tuan](https://github.com/chaoi111) (original development) and interns of i5k workspace.
## Requirement
* python2.7
* [EDirect](https://www.ncbi.nlm.nih.gov/books/NBK179288/#chapter6.Getting_Started)
## Background
* We need to submit the Official Gene Sets generated at the i5k Workspace to NCBI/GenBank. To do so, we need to 1) clean up of the original gff, and 2) run NCBI's new table2asn program. 
* [More information on NCBI's requirements for the gff3 file](https://www.ncbi.nlm.nih.gov/sites/genbank/genomes_gff/)
## Usage
    python  ogs_v5.py  assembly_ID  gff3_file
    example: ~/anaconda2/bin/python ogs_v5.py GCA_000696205.1 clec_OGS_v1_2_with_pep_CDS.2.gff3


## Current Functions (Function marked with strikethrough has been removed)
#### 2. Remodel_pseudogenes()
    - change parent feature type from 'pseudogene' to 'gene'     
#### 3. check_locus()
    - scan gff file to see whether some, all, or no gene features have the attribute 'locus_tag' in column 9 (type gene)   
#### 3. locus_tag_util(Biosample) 
    - identify the locus tag prefix associated with the genome project's BioProject   
#### 4. ID_attribute()  
    - create transcript_id attribute for mRNA, format: gnl|[db]|[ID]      
#### 4. cds_id_attribute()   
    - create protein_id attribute for CDS, format: gnl|[db]|[ID] 
#### 5. ~~remove_Note()~~
    - replace Note attribute value with "Additional notes are available on the i5k Workspace gene page."
#### 6. ~~name_Attributes_gene()~~ 
    - copy mRNA feature 'Name' attribute as a 'product' attribute
#### 7.~~name_Attribute_mRNA_pseudo()~~ 
    - replace Name attribute with product attribute
#### 8. ~~name_Attribute_trncRNA()~~  
    - replace Name attribute with product attribute
#### 9. Scan_for_genome_alterations() 
    -  output summary line including number of features found of each type
#### 10. ~~description_attribute()~~  
    -  scan for description attributes in features other than gene features
#### 11. ~~symbol_attribute()~~ 
    - scan for symbol attributes in features other than gene features
#### 12. Dbxref()  
    - gene features: Replace all existing Dbxrefs with this Dbxref: I5KNAL:[ID], where [ID] is the ID attribute of the gene. 
    - all other features: remove all Dbxrefs.
