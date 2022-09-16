import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-sam', required=True, help='sam file')
args = vars(parser.parse_args())
sam_file = args['sam']

# sam_file = '/Users/songweizhi/Desktop/subset.sam'
# sam_file = 'combined_pass_rd2_trimmed50bp_first_1000bp.sam'

mapped_num = 0
unmapped_num = 0
for each in open(sam_file):
    if not each.startswith('@'):
        each_split = each.strip().split('\t')
        ref_id = each_split[2]
        if ref_id == '*':
            unmapped_num += 1
        else:
            mapped_num += 1

print('%s\t%s/%s\t%s' % (sam_file, mapped_num, (mapped_num + unmapped_num), (mapped_num*100/(mapped_num + unmapped_num))))
