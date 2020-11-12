from celescope.__init__ import __CONDA__
from celescope.rna.__init__ import __STEPS__, __ASSAY__
from celescope.tools.Multi import Multi


class Multi_rna(Multi):
    pass

def main():
    multi = Multi_rna(__ASSAY__, __STEPS__, __CONDA__)
    multi.run()

if __name__ == '__main__':
    main()
