from witch_msa.configs import Configs
from witch_msa.helpers.alignment_tools import ExtendedAlignment
import gzip

'''
Callback for a query alignment result
*args should have three fields: _query, _index, checkpoint_path
'''
def callback_queryAlignment(success, ignored, retry, i_retry, 
        query, index, taxon_name, checkpoint_path):
    if (not query) and i_retry > 0:
        retry.append(index)
    else:
        # failed job indicated in the <witch-ng> or <default> pipelines,
        # should be ignored in the output
        if (not query) or len(query) == 0:
            ignored.append(taxon_name)
        else:
            # write to checkpoint_path
            if (not isinstance(query, ExtendedAlignment)) or (len(query) != 1):
                return
            seq = query[taxon_name]
            line = '{}\t{}\n'.format(taxon_name, seq)
            encoded = line.encode('utf-8')
            with gzip.open(checkpoint_path, 'ab') as f:
                f.write(encoded)

            # update success
            success.append(query)
