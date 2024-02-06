import os

class Resources:
    """Gets URLs and file names of fasta and GTF files for a given genome and build
    """
    
    # create genome directory
    os.makedirs("resources/", exist_ok=True)
    
    def __init__(self, genome, build):
        self.genome = genome
        self.build = build
                
        # base URLs
        base_url = f"https://ftp.ensembl.org/pub/release-{build}/"
                
        if "hg" in genome:
            if genome == "hg19":
                name = "GRCh37"
                                
                self.blacklist_url = "https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg19-blacklist.v2.bed.gz"
            elif genome == "hg38":
                name = "GRCh38"
                                                
                self.blacklist_url = "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz"
                
            # create URLs for genome files
            self.fasta_url = f"{base_url}fasta/homo_sapiens/dna/Homo_sapiens.{name}.dna.primary_assembly.fa.gz"
            self.mt_fasta_url = f"{base_url}fasta/homo_sapiens/dna/Homo_sapiens.{name}.dna.chromosome.MT.fa.gz"
            self.gtf_url = f"{base_url}gtf/homo_sapiens/Homo_sapiens.{name}.{build}.gtf.gz"
            
        elif "mm" in genome:
            if genome == "mm38":
                name = "GRCm38"
                
            elif genome == "mm39":
                name = "GCRm39"
                               
            # create URLs for genome files
            self.fasta_url = f"{base_url}fasta/mus_musculus/dna/Mus_musculus.{name}.dna.primary_assembly.fa.gz"
            self.mt_fasta_url = f"{base_url}fasta/mus_musculus/dna/Mus_musculus.{name}.dna.chromosome.MT.fa.gz"
            self.gtf_url = f"{base_url}gtf/mus_musculus/Mus_musculus.{name}.{build}.gtf.gz"
        
        elif "dm" in genome:
            if genome == "dm6":
                name = "BDGP6.46"

            self.fasta_url = f"{base_url}fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.{name}.dna.toplevel.fa.gz"
            self.mt_fasta_url = f"{base_url}fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.{name}.dna.MT.fa.gz"
            self.gtf_url = f"{base_url}gtf/drosophila_melanogaster/Drosophila_melanogaster.{name}.{build}.gtf.gz" 
        
        # downloaded unzipped file names
        #self.fasta = self._file_from_url(self.fasta_url)
        #self.gtf = self._file_from_url(self.gtf_url)
        
        elif genome == "MG1655":
            # Spike-in genome
            self.fasta_url = f"https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-{build}/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/dna/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.dna.toplevel.fa.gz"
            self.gtf_url = "" # placeholder
            self.blacklist_url = "" # placeholder
            self.nomt_fasta_url = "" # placeholder
            self.mt_fasta_url = "" # placeholder
        
        # downloaded unzipped file names
        self.fasta = self._file_from_url(self.fasta_url)
        self.gtf = self._file_from_url(self.gtf_url)
        self.blacklist = self._file_from_url(self.blacklist_url)
        self.mt_fasta = self._file_from_url(self.mt_fasta_url)
        self.nomt_fasta = self.fasta.replace(".fa", ".no_mt.fa")
        
    def _file_from_url(self, url):
        """Returns file path for unzipped downloaded file
        """
        
        return f"resources/{os.path.basename(url).replace('.gz','')}"
      
        
            
            
    