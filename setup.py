from distutils.core import setup, Extension
from distutils.command.build_ext import build_ext as DistutilsBuild
import subprocess


class MyBuild(DistutilsBuild):
    def install_libs(self):
       if(subprocess.call("cd libs && make",shell=True)):
           raise RuntimeError,'Make of seq-align failed'

    def run(self):
       self.install_libs()
       DistutilsBuild.run(self)

module1 = Extension('seqalign', sources=['seqalign.c'],
                    include_dirs = ['libs/seq-align/src'],
                    libraries =['align'],
                    library_dirs = ['libs/seq-align/src'])
setup(
    name="seqalign",
    version="0.1",
    description= 'Wrapper for sequence alignment',
    py_modules = ["seqal"],
    ext_modules = [module1],
    cmdclass={'build_ext':MyBuild})

