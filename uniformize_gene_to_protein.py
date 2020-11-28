import re
MGG_rgx=re.compile('^MGG_')
def get_id_protein(old_id):
    if MGG_rgx.match(old_id) is not None:
        return old_id
    else:
        parts=old_id.split('_')
        if len(parts)==5:
            return '{}_{}'.format(parts[0],parts[3])
        else:
            assert len(parts)==4
            return '{}_{}'.format(parts[0],parts[2])
        # elif len(parts)==4:
        #     return '{}_{}'.format(parts[0],parts[2])
        # else:
        #     return old_id
def get_id_gene (old_id):
    if MGG_rgx.match(old_id) is not None:
        return old_id+'T0'
    else:
        parts=old_id.split('_')
        if len(parts)==5:
            return '{}_{}'.format(parts[0],parts[3])
        else:
            assert len(parts)==4
            return '{}_{}'.format(parts[0],parts[2])