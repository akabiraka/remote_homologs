import sys
sys.path.append("../remote_homologs")

f = open("data/generated/cdhit_clusters/scop_sf_represeq_at_90_seq_identity.txt", mode="r")
out = open("data/generated/cdhit_clusters/scop_sf_represeq_lib_latest.35clust", mode="w")

flag=False
for i, line in enumerate(f.readlines()):
    # print(i, line)
    
    if ">Cluster" not in line:
        line_items = line.split(" ")
        cluster.append(str(line_items[1][1:-3]))
        print(cluster)

    elif ">Cluster" in line:
        if flag: 
            out.writelines(" ".join(cluster))
            out.write("\n")
        
        cluster = []
        flag = True
        # print(" ".join(cluster))
        

    
    # if i==5000: break
