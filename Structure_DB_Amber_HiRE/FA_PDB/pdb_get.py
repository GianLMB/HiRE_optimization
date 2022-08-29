#! /usr/bin/env python

import urllib, re, os, sys, getopt

def usage():
  print "\nusage: pdb_get [options] <code> "
  print "    where [options] could be:"
  print "       -p to retrieve PDB format (default)"
  print "       -c to retrieve mmCIF format"
  print "       -s to retrieve structure factors along with the PDB format coordinates"
  print "       and <code> is the 4-character PDB entry code"

def get_options():
  pdb = 1
  mmCIF = 0
  struct_fact = 0
  try:
    opts,args = getopt.getopt(sys.argv[1:],'hpcs')
  except:
    print 'Unrecognized Option: ', sys.argv[1:]
    usage()
    return pdb, mmCIF, struct_fact

  for o,a in opts:
    if o == '-h':
      usage()
      sys.exit(0)
    elif o == '-p':
      pdb = 1
    elif o == '-c':
      pdb = 0
      mmCIF = 1
    elif o == '-s':
      struct_fact = 1

  return pdb, mmCIF, struct_fact, args

def main():
  (pdb, mmCIF, struct_fact, args) = get_options()

  for code in args:
    code = code.lower()

    if ( pdb == 1):
      print "\nDownloading %s.pdb.gz .........." % (code),
      url = 'ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb%s.ent.gz' % code
      filename = code + '.pdb.gz'
      try:
        urllib.urlretrieve(url, filename)
        print "Uncompressing %s.pdb.gz" % code
        os.system("gunzip %s.pdb.gz" % code)
      except:
        print "Error retrieving %s" % url
        
    elif (mmCIF == 1):
      print "\nDownloading %s.cif.gz .........." % (code),
      url = 'ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/mmCIF/%s.cif.gz' % code
      filename = code + '.cif.gz'
      try:
        urllib.urlretrieve(url, filename)
        print "Uncompressing %s.cif.gz" % code
        os.system("gunzip %s.cif.gz" % code)
      except:
        print "Error retrieving %s" % url

    if ( struct_fact == 1 ):
      print "\nDownloading r%ssf.ent.gz .........." % (code),
      url = 'ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/structure_factors/r%ssf.ent.gz' % code
      filename = 'r' + code + 'sf.ent.gz'
      try:
        urllib.urlretrieve(url, filename)
        print "Uncompressing r%ssf.ent.gz" % code
        os.system("gunzip r%ssf.ent.gz" % code)
      except:
        print "Error retrieving %s" % url

if __name__ == '__main__':
  main()
