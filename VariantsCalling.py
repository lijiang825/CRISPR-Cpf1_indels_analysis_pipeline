import math,sys

def samParser(sam_file):
    '''Parse a sam file; discard headers;
        only keep sequence(0), flag(1), reference(2) and CIGAR(3) 
    '''
    parsed = []
    for i in sam_file:
        temp = i.strip()
        if temp[0] != "@":
            temp2 = temp.split("\t")
            parsed.append((temp2[9],temp2[1],temp2[2],temp2[5]))
    return parsed

def cigarParser(cigar_string):
    '''Parse a cigar string into individual deletion/insertion: [(mut,size),...]
    '''
    parsed = []
    nt_size = ""
    for i in cigar_string:
        if i == "M" or i == "I" or i == "D" or i == "S" or i == "H":
           parsed.append((i,int(nt_size)))
           nt_size = ""
        else:
            nt_size += i
    return parsed

def indelInGuide(M,D,dToGuide,ref_orientation):
    '''Compute nt of deletions in the guide sequence + PAM
        M,D: matched and deleted bases from parsed SAM string; M always comes before D
        dToGuide: distance from start of sequence to 3' of guide RNA
        dToGuide(GFP_R49) = 48
        dToGuide(plasmid) = 120
        ref_orientation: 3'-5' = -1; 5'-3' = 1 
    '''
    if ref_orientation == -1:
        minRange, maxRange = dToGuide, dToGuide + 24
    else:
        minRange, maxRange = dToGuide-23, dToGuide+1
    indelInGuide = 0
    for i in range(M+1,M+D+1):
        if i in range(minRange,maxRange):
            indelInGuide += 1
    return indelInGuide

def cigarProcess(cigar_parsed,dToGuide,ref_orientation):
    ''' Extract indel identity, size and middle position to 3' of guide RNA from cigar mutation and size list
        Store in a dict and return the dict
        dToGuide(GFP_R49) = 48
        dToGuide(plasmid) = 120
    '''
    ## Edited = 0 (not edited); 1 (edited) 2 (discarded due to wrong cigar string)
    result = {"edited":0, "indelSize":[], "indelPos":[], "indelGuide":[], "truePos":[]} 
    
    ## If cigar_parsed is empty, discard this cigar
    if len(cigar_parsed) == 0:
        result["edited"] = 2
        return result
    
    ## If cigar has "S", "H" or first element is not "M", discard this cigar
    ## Specifically for BWA
    if cigar_parsed[0][0] != "M":
        result["edited"] = 2
        return result
    for i in xrange(len(cigar_parsed)):
        if cigar_parsed[i][0] == "S" or cigar_parsed[i][0] == "H":
            result["edited"] = 2
            return result

    ## If mut[-1] == "D" and size[-1] >= 100; remove the last elements
    ## Specifically for pairwise aligner 
    if cigar_parsed[-1][0] == "D" and cigar_parsed[-1][1] >= 100: 
        cigar_parser2 = cigar_parsed[:-1]
    else:
        cigar_parsed2 = cigar_parsed
        
    ## Scan the mut/size list until the last "M" with size > =20; make new mut/size lists (3) with that subset
    ## Specifically for pairwise aligner 
    lastM = 0
    for i in xrange(len(cigar_parsed2)):
        if cigar_parsed2[i][0] == "M" and cigar_parsed2[i][1] >= 20:
            lastM = i
    cigar_parsed3 = cigar_parsed2[:lastM+1]
    
    ## Compute and store the size, true position, relative position to 3' of guide RNA,
    ## and size of guide RNA removed by edits into the dict 
    cigar_length = len(cigar_parsed3)
    ## Return WT without any edit
    if cigar_length == 1:
        result["edited"] = 0
        return result
    else:
        result["edited"] = 1
    ## pointer for keeping track of true position 
    pointer = 0 
    for i in xrange(cigar_length):
        current_type, current_size = cigar_parsed3[i][0], cigar_parsed3[i][1]
        ## Skip match, only advance pointer
        if current_type == "M":
            pointer += current_size
            continue
        ## Store true position
        result["truePos"].append(pointer)
        ## Compute indel size
        if current_type == "D":
            result["indelSize"].append(-1 * current_size)
        if current_type == "I":
            result["indelSize"].append(current_size)
        ## Compute relative position to 3' of guide RNA
        result["indelPos"].append(ref_orientation * (pointer + int(math.ceil(0.5 * current_size))-dToGuide))
        ## Compute deletion in guide RNA
        if current_type == "D":
            result["indelGuide"].append(indelInGuide(pointer,current_size,dToGuide,ref_orientation))
        else:
            result["indelGuide"].append(0)
        ## advance pointer when type=="D" 
        if current_type == "D":
            pointer += current_size
    return result

def buildIndel(indelSize, truePos):
    if indelSize > 0:
        result = str("I") + str(truePos)
    else:
        result = str("D") + str(truePos)
    return result

def main(argv):
    ## Make connection for input and output files
    inputPath = argv[0]
    inputSAM = open(inputPath)

    ## Parse SAM file and initiate result and meta result dict
    parsedSam= samParser(inputSAM)
    par_dict = {"eGFP":(48,-1),"Cpf1_DR_gRNA_vector_rev":(120,1)}
    result = {"eGFP":[], "Cpf1_DR_gRNA_vector_rev":[]}
    ## [nTotal,nMut,nDel,nIns]
    result_sum = {"eGFP":[0]*4,"Cpf1_DR_gRNA_vector_rev":[0]*4}
    indel_cigar_set = {"eGFP": set(), "Cpf1_DR_gRNA_vector_rev": set()}
  
    
    for item in parsedSam:
        ## Check flag
        if item[1] != "0":
            continue
        ## Parse cigar string
        ref = item[2]
        cigar_parsed = cigarParser(item[3])
        cigar_dict = cigarProcess(cigar_parsed, par_dict[ref][0],par_dict[ref][1])
        ## Check 'Edited' flag 
        if cigar_dict["edited"] == 2:
            continue
        if cigar_dict["edited"] == 0:
            result_sum[ref][0] += 1
        if cigar_dict["edited"] == 1:
            result_sum[ref][0] += 1
            check = 0
            size = len(cigar_dict["indelSize"])
            if len(cigar_dict["indelPos"]) != size or len(cigar_dict["indelGuide"]) != size or len(cigar_dict["truePos"]) != size:
                raise ValueError("All items in the cigar_dict should have identical lengths")
            for i in xrange(size):
                if abs(cigar_dict["indelSize"][i]) < 3 and abs(cigar_dict["indelPos"][i]) > 25:
                    continue
                else:
                    check = 1                                                         
                    result[ref].append([item[3], buildIndel(cigar_dict["indelSize"][i], cigar_dict["truePos"][i]),
                                       cigar_dict["indelSize"][i], cigar_dict["indelPos"][i], cigar_dict["indelGuide"][i]])
                    if cigar_dict["indelSize"][i] < 0:
                        result_sum[ref][2] += 1
                    else:
                        result_sum[ref][3] += 1
            result_sum[ref][1] += check
            if check == 1:
                indel_cigar_set[ref].add(item[3])

    sample = inputPath.split('/')[-1]
    print "====================================="
    print "Print results for " + sample
    print "************************************"
    for key in result_sum:
        print "Results for " + key
        print "Total number of reads: " + str(result_sum[key][0]) 
        print "Total number of reads with indels: " + str(result_sum[key][1])
        result_sum[key].append(float(result_sum[key][1])/result_sum[key][0])
        print "The frequency of indel is ", round(result_sum[key][4],4)
        print "Total number of reads with deletions: " + str(result_sum[key][2]) 
        print "Total number of reads with insertions: " + str(result_sum[key][3])
        result_sum[key].append(len(indel_cigar_set[key]))
        print "Total number of unique indels: " + str(result_sum[key][5])
        print "************************************"
    print "====================================="

    for ref in result_sum:
        outputPath = inputPath+"."+ref+".csv"
        outputCSV = open(outputPath,"w")
        outputCSV.write("cigar,indel,size,position,delGuide\n")
        for entry in result[ref]:
            temp = ""
            for column in entry:
                temp = temp + str(column) + ","
            outputCSV.write(temp.rstrip(",") + "\n")

        outputPath = inputPath+"."+ref+".meta"
        outputCSV = open(outputPath,"w")
        outputCSV.write("sumReads,sumMutations,sumDeletions,sumInsertions,mutationFreq,uniqueIndel\n")
        temp = ""
        for entry in result_sum[ref]:
            temp = temp + str(entry) + ","
        outputCSV.write(temp.rstrip(",") + "\n")

    
    inputSAM.close()
    outputCSV.close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise ValueError("Usage: input, sam file of alignment from BWA")
    main(sys.argv[1:])


