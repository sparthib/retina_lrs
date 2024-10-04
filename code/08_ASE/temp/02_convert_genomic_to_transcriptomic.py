import gzip
import sys
from copy import deepcopy

sys.path.insert(0, './')
from ref_lib.GTF import GTFfile, GTFEntry, get_gtf_contents
from ref_lib.Fasta import FastaEntry, FastaFile

from collections import defaultdict, OrderedDict


reverse_complement_dict = { "A" : "T",
                            "C" : "G",
                            "G" : "C",
                            "T" : "A",
                            "N" : "N"}

VCF_FIELDS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",  "INFO",    "FORMAT",  "CAST_EiJ"]

class VcfEntry:

    def __init__(self , vcf_line_contents ):
        assert len(vcf_line_contents) >= len(VCF_FIELDS)
        
        self.fields = { VCF_FIELDS[i] : vcf_line_contents[i] for i in range( len(VCF_FIELDS) ) }
        
        
    def __str__(self ):
        """
        This needs to be rewritten 
        """
        return "\t".join( [self.fields[f] for f in VCF_FIELDS] )


############################################################################
    
class VcfFile:
    '''
    This is a reader for 
    '''
    
    def __init__(self , file):
        myopen = open
        if file.endswith(".gz"):
            myopen = gzip.open

        if(file):
            self.f = myopen(file , "rt")
        else:
            self.f = stdin

    #####################################################

    def __enter__(self):
        return self

    #####################################################

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    ######################################################

    def __getitem__(self, index):
        line = self.f.readline().strip()
        
        while line.startswith("#"):
            line = self.f.readline().strip()
        
        if line == "":
            raise IndexError
        #line_contents = line.split("\t")
        line_contents = line.split()
        if len(line_contents) < 9:
            raise IndexError
        return VcfEntry(line_contents)
                
    #########################################################

    def __del__(self):
        self.f.close()
        
vcf_file = "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/01_longshot_vcfs/H9-FT_1/H9-FT_1.vcf"

my_vcf = VcfFile(vcf_file)


gtf_file = "/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf"
gtf_all = get_gtf_contents(gtf_file)

fasta_file   = "/dcs04/hicks/data/sparthib/references/transcriptome/GENCODE/gencode.v44.transcripts.fa" 
fasta_reader = FastaFile(fasta_file)


fasta_entries = [a for a in fasta_reader]

fasta_dict = OrderedDict()

for f in fasta_entries:
    fasta_dict[f.header] = f.sequence
    
# We find the fasta entries of the transcript files

fasta_gtf_entries = {}
transcript_to_header = {}

for f in fasta_entries:
    contents                    = f.header.split("|")
    gene_id                     = contents[1].split(".")[0]
    transcript_id               = contents[0].split(".")[0]
    fasta_gtf_entries[f.header] = gtf_all[gene_id][transcript_id]
    
    transcript_to_header[transcript_id] = f.header
    
    
# Now we convert + strand coordinates
conversion_dict = dict()

for fasta_header, t_dict  in fasta_gtf_entries.items():
    
    if not conversion_dict.get(t_dict["chr"]):
        conversion_dict[t_dict["chr"]] = defaultdict(list)
        
    transcript_id = fasta_header.split("|")[0].split(".")[0]    
        
    t_position = 1
        
    if t_dict["strand"] == "+":
        
        temp_dict     = defaultdict(list) 
        
        for e in t_dict["exons"]:
            for i in range(e[0], e[1] + 1):
                conversion_dict[t_dict["chr"]][i].append( (transcript_id, t_position, "+" ) )
                t_position += 1
    else:
        for e in t_dict["exons"]:
            for i in range(e[1], e[0] - 1, -1):
                conversion_dict[t_dict["chr"]][i].append( (transcript_id, t_position, "-" ) )
                t_position += 1

for pos, mylist in conversion_dict['chr4'].items():
    if (len(mylist) > 1):
        print("{} : {}".format(pos, mylist) )
        
counter = 0

f_handle = gzip.open("transcriptomic_variants.vcf.gz", "wt")

for v in my_vcf:
    if conversion_dict[v.fields['CHROM']].get( int(v.fields['POS']) ):
        counter +=1
        
        for entry in  conversion_dict[v.fields['CHROM']][ int(v.fields['POS']) ]: 
            this_entry = deepcopy(v)
            this_entry.fields['CHROM'] = transcript_to_header[entry[0]]
            this_entry.fields['POS']   = str(entry[1])
            
            if entry[2] == "-":
                this_entry.fields['REF'] =  "".join(map( lambda x: reverse_complement_dict[x], this_entry.fields['REF']))
                this_entry.fields['ALT'] =  "".join(map( lambda x: reverse_complement_dict[x], this_entry.fields['ALT']))
                
            print(this_entry, file = f_handle)
            
 
f_handle.close()
        
print("Variants found {}".format(counter))

## Sanity_cehck
## Are reference nucleotides in the Vcf file same as the nucleotides in the fasta files?

produced_vcf_file = "transcriptomic_variants.vcf.gz"
#produced_vcf_file = "mock.vcf.gz"

produced_vcf = VcfFile(produced_vcf_file)


# This should NOT print anything on the produced_vcf_file
# We put one fault record on the mock file to test our test :)

for v in produced_vcf:
    observed_nuc = v.fields["REF"]
    
    # Note that vcf is 1-based and lists are 0-based
    # so we need -1 in the index
    expected_nuc = fasta_dict[v.fields["CHROM"]][ int(v.fields["POS"])-1]
    
    if observed_nuc != expected_nuc:
        print("In", v.fields["CHROM"])
        print("At position ", v.fields["POS"])
        print("{} != {}".format(observed_nuc, expected_nuc) )
        

## Now we mask our Transcriptomic fasta file based on the variants on VCF

produced_vcf = VcfFile(produced_vcf_file)

masked_fasta_dict = deepcopy(fasta_dict)

for v in produced_vcf:
    #print(v.fields["POS"])
    #print(v.fields["CHROM"])
    #print(masked_fasta_dict[v.fields["CHROM"]])
    masked_sequence = masked_fasta_dict[v.fields["CHROM"]][ :(int(v.fields["POS"]) - 1)  ] + \
                      "N" + \
                      masked_fasta_dict[v.fields["CHROM"]][ int(v.fields["POS"]): ]
    masked_fasta_dict[v.fields["CHROM"]] = masked_sequence
    


with gzip.open("variant_masked_transcriptome.fa.gz", "wt") as out_stream:
    for header, sequence in masked_fasta_dict.items():
        this_entry = FastaEntry(header = header, sequence = sequence)
        print(this_entry, file = out_stream)
        

# Let's make sure that all masked positions are N's

masked_fasta_file   = "variant_masked_mouse_transcriptome.fa.gz"
masked_fasta_reader  = FastaFile(masked_fasta_file)

read_fasta_dict = {}

for f in masked_fasta_reader:
    read_fasta_dict[f.header] = f.sequence

produced_vcf_file = "transcriptomic_variants.vcf.gz"
#produced_vcf_file = "mock.vcf.gz"

produced_vcf = VcfFile(produced_vcf_file)

for v in produced_vcf:
    
    # Note that vcf is 1-based and lists are 0-based
    # so we need -1 in the index
    assert read_fasta_dict[v.fields["CHROM"]][ int(v.fields["POS"])-1] == "N"


