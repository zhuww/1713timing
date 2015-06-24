from datatools.tempo import *

tf = TOAfile('1713.Feb.T2.tim')
print tf.groups.keys()

ntf = tf.subgroup(groups=['M3A-L', 'M3B-L', 'ABPP-L', 'ABPP-S', 'M4-L', 'M4-S'])
foo = open('1713.sns.tim', 'w')
foo.write(ntf.tempo2fmt())
foo.close()
