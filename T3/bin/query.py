import cPickle as pickle
from os import getenv


def query():
      r = []
      try:
            l = pickle.load(open(getenv('SUBMIT_WORKDIR')+'/submission.pkl'))
            s = l[-1]
            r.append( 'ClusterID '+str(s.cluster_id) )
            statii = s.query_status()
            r.append( 'Job summary:' )
            for k,v in statii.iteritems():
                  r.append( '\t %10s : %5i'%(k,len(v)) )
            return r
      except IOError:
            r.append( 'No job submitted yet!' )
            return r

if __name__ == '__main__':
      print '\n'.join(query())
