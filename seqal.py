from seqalign import *

#copied from nwalign module
def read_matrix(path):
    fh = open(path)
    headers = None
    while headers is None:
        line = fh.readline().strip()
        if line[0] == '#': continue
        headers = [x for x in line.split(' ') if x]

    line = fh.readline()
    result = []
    while line:
        h1 = line[0]
        line = [int(x) for x in line[1:-1].split(' ') if x]
        result.extend([((h1,h),l) for h,l in zip(headers,line)])
        line = fh.readline()
    return dict(result)

def alignment_to_mapping(s1, s2):#{{{
    mapping = []
    cur_value = None 
    for c1, c2 in zip(s1,s2):
        x1 = c1 != '-'
        x2 = c2 != '-'
        if x1 and x2:
            if not isinstance(cur_value,int):
                mapping.append(cur_value)
                cur_value = 0
            cur_value += 1
        elif x1:
            if not isinstance(cur_value, list) or cur_value[0] == 0:
                mapping.append(cur_value)
                cur_value = [0,0]
            cur_value[0] += 1
        elif x2:
            if not isinstance(cur_value, list) or cur_value[1] == 0:
                mapping.append(cur_value)
                cur_value = [0,0]
            cur_value[1] += 1
        else:
            raise RuntimeError, 'Double gap in %s and %s' % (s1,s2)
    mapping.append(cur_value)
    mapping = mapping[1:] #remove none
    pos = 0
    while pos < (len(mapping) - 1):
        if isinstance(mapping[pos],list) and isinstance(mapping[pos+1], list):
            x = mapping[pos]
            x[0] += mapping[pos+1][0]
            x[1] += mapping[pos+1][1]
            del mapping[pos+1]
        else:            
            pos += 1


    mapping = tuple([tuple(m) if isinstance(m,list) else m for m in mapping])
    return mapping#}}}
