# caspr_search
Dependencies:

  vmatch2: http://vmatch.de/distributions/vmatch-2.3.0-Linux_x86_64-64bit.tar.gz
  
  emboss: sudo apt-get install emboss
  
  prodigal: https://github.com/hyattpd/Prodigal
  
  hmmer: http://eddylab.org/software/hmmer/hmmer-3.3.tar.gz
  
  biopython: pip3 install biopython

usage:
python3 -i test.fasta -o out_dir --meta
