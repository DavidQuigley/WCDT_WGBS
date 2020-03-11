samples = ['DTB-003-BL','DTB-005-BL','DTB-008-BL','DTB-009-BL','DTB-011-BL','DTB-018-BL','DTB-019-PRO','DTB-020-BL', \
 'DTB-021-BL','DTB-022-BL','DTB-023-BL','DTB-024-PRO','DTB-034-BL','DTB-035-BL','DTB-036-BL','DTB-037-BL', \
 'DTB-040-BL','DTB-042-BL','DTB-053-BL','DTB-055-PRO','DTB-059-BL','DTB-060-BL','DTB-061-BL','DTB-063-BL', \
 'DTB-064-BL','DTB-067-PRO','DTB-069-BL','DTB-071-BL','DTB-074-BL','DTB-077-PRO','DTB-080-BL','DTB-083-BL', \
 'DTB-085-BL','DTB-089-BL','DTB-090-PRO','DTB-091-BL','DTB-092-BL','DTB-094-BL','DTB-097-PRO','DTB-098-PRO2', \
 'DTB-100-BL','DTB-101-BL','DTB-102-PRO','DTB-104-BL','DTB-111-PRO','DTB-112-BL','DTB-119-PRO','DTB-121-BL', \
 'DTB-124-BL','DTB-126-BL','DTB-127-PRO','DTB-128-BL','DTB-129-BL','DTB-132-BL','DTB-135-PRO','DTB-137-PRO', \
 'DTB-138-BL','DTB-140-BL','DTB-141-BL','DTB-143-BL','DTB-146-BL','DTB-149-BL','DTB-151-BL','DTB-156-BL', \
 'DTB-159-BL','DTB-165-PRO','DTB-167-PRO','DTB-170-BL','DTB-172-BL','DTB-173-BL','DTB-175-BL','DTB-176-BL', \
 'DTB-183-BL','DTB-186-BL','DTB-187-BL','DTB-188-BL','DTB-190-BL','DTB-194-PRO','DTB-201-PRO','DTB-202-BL', \
 'DTB-204-BL','DTB-205-BL','DTB-206-BL','DTB-210-BL','DTB-213-BL','DTB-214-BL','DTB-216-PRO','DTB-220-BL', \
 'DTB-222-BL','DTB-223-BL','DTB-232-PRO','DTB-234-BL','DTB-251-BL','DTB-252-BL','DTB-255-BL','DTB-258-BL', \
 'DTB-260-BL','DTB-261-BL','DTB-265-PRO','DTB-266-BL']

chroms = ['NC_000001.11','NC_000002.12','NC_000003.12','NC_000004.12','NC_000005.10','NC_000006.12', \
'NC_000007.14','NC_000008.11','NC_000009.12','NC_000010.11','NC_000011.10','NC_000012.12','NC_000013.11', \
'NC_000014.9','NC_000015.10','NC_000016.10','NC_000017.11','NC_000018.10','NC_000019.10','NC_000020.11', \
'NC_000021.9','NC_000022.11','NC_000023.11','NC_000024.10']

fn_chrom = "/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/metadata/coverage/cpg_coverage_per_chromosome_both_strands_test.txt"
fn_samples = "/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/metadata/coverage/cpg_depth_by_sample_both_strands_test.txt"
fo_chrom = open(fn_chrom, 'w')
fo_samples = open(fn_samples, 'w')
dir = "/data1/datasets_1/human_prostate_WCDT/wgbs_solid/processed/bismark/cpg_reports/tumors"

fo_chrom.write( "chrom\t" + '\t'.join( ["s" + x for x in [str(x) for x in range(1,101)]] ) + '\n' )
fo_samples.write( "sample_id\tchrom\tn_0\tn_10\tn_20\tn_30\tn_40\n")
for chrom in chroms:
    cpg_cover2sample = {}
    sample2coverage = {}
    for sample in samples:
        n_10 = 0
        n_20 = 0
        n_30 = 0
        n_40 = 0
        n_0 = 0
        fn = dir + "/" + sample + ".merged.CX_report.txt.gz." + chrom + ".methyl"
        print( chrom + " " + sample )
        f = open(fn)
        pos_prev = -100
        for line in f:
            a = line.rstrip('\r\n').split('\t')
            depth = int(a[2])
            pos = int( a[0] )
            if pos == pos_prev+1:
                n_0 += 1
                depth_total = depth + depth_prev
                if depth_total >= 10:
                    n_10 += 1
                    try:
                        cpg_cover2sample[ pos_prev ] = cpg_cover2sample[ pos_prev ] + 1
                    except KeyError:
                        cpg_cover2sample[ pos_prev ] = 1
            
                if depth_total >= 20:
                    n_20 += 1
            
                if depth_total >= 30:
                    n_30 += 1
            
                if depth_total >= 40:
                    n_40 += 1
            
            pos_prev = pos
            depth_prev = depth
        
        sample2coverage[sample] = [n_0, n_10, n_20, n_30, n_40]
    
    covers = [0] * 100
    for k in cpg_cover2sample:
        covers[ cpg_cover2sample[k] - 1 ] += 1
    
    fo_chrom.write( chrom + '\t' + '\t'.join( [str(x) for x in covers] ) + '\n' )
    for sample in samples:
        fo_samples.write( sample + '\t' + chrom + '\t' + '\t'.join( [str(x) for x in sample2coverage[sample] ] ) + '\n' )

fo_samples.close()
fo_chrom.close()
