# Ideas for the depth smoothing and binning



# 0. Smooth outliers by first replacing values with a (small) window mean. Maybe 5mers? - smoothwin=X
# 1. Set a number of fixed points and a number of flexible points for k means clustering - binnum=X fixedbin=LIST
# 2. Perform k means clustering to bin sites - bincycle=NUM 
# 3. Compress vector into regions (bed file)
# 4. Smooth out the small regions?



# Set a vector of L-1 as the boundaries, and make TRUE where dep[i] != dep[i+1]
# 
# » Still need to extract lengths and smooth the small ones out. (Default = 100 bp region?)
# 
# » Want to look for the transitions from diploid to haploid and haploid to lowcov. Will want to identify collapsed regions but possibly ignore them if using for a new Diploidocus mode
# 
# » Pull out the non-diploid regions and map then look for non-self hits that are also sub-diploid? (Check out the purge_dups settings.)
# 
# » Also need to identify and deal with gaps. Gap-adjacent regions can be trimmed more aggressively that internal. Mark gaps as -1 depth rather than 0 and make own category. Incorporate them into diploid or haploid regions? Or split regions on gaps? (Maybe split to start with.)
# 
# Q. Should we allow overlapping regions, keep small unassigned regions, or arbitrarily split mid-region if a small one is flanked by two different ones? (Maybe merge based on the distance, rather than filter based on the length? Would take care of flip-flopping. Give diploid priority to avoid issues of alteranting sequences.



# Name ideas - want something fishy
# Region Depth Bin Classification Genome Assembly
# GARDian
# dogfish? Seems to be part of nanopore tools. 
# DOGSHARC - Depth of Genome Sequencing Haploid Assembly Region Classifier

# GenomeARC = old idea
# DOGSARC