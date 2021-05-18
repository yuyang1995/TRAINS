import re
import os


def parse(fname, parset):
    regex = r'^(\S+)\s*(\S+)'
    if os.path.exists(fname):
        with open(fname) as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if not line.startswith('%'):
                    mo = re.search(regex, line)
                    if mo:
                        parset[mo.group(1)] = mo.group(2)


parset = {}
parse('../param/param', parset)
parse('../' + parset['ParamPT'], parset)
parse('../data/DNEST_OPTIONS', parset)
