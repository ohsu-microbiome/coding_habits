import unittest
import os
import subprocess
from datetime import datetime
from Bio import AlignIO
from Bio import SeqIO
import shutil


class test_ambitious(unittest.TestCase):
    """
        subprocess.run() is new to me. Call the function and any 
        options as a list of arguments. 

        I had a lot of problems at first but figured it out.
            
        The paths to the different files were confusing. 
        Settled on using the cwd="where/is/foo.bar" to change the working directory 
        to the right /src_files dir, then point to the testing file by relative path.
            
        BioPython got confused when an argument in the subprocess list was like this:
            
            subprocess.run(["python", "extract_demo.py", "-i bad_clustal_format.aln"]
                                                                  ^
            
        There's something going on in the Windows environment that makes the OS 
        look for a file named " foo.bar", a file with a space as the first character.
        Changing that to split up the flag from the path fixed that:
            
            subprocess.run(["python", "extract_demo.py", "-i", "bad_clustal_format.aln"]
                                                               ^
                    
        Now you know.
    """
    def test_script_exists(self):
        """
            check the script exists
        """
        get_files=os.listdir("./")
        self.assertIn("extract_demo.py", get_files)
        
    def test_find_resources(self):
                
        get_files=os.listdir("./")
        self.assertIn("good_clustal_format.aln", get_files)
        self.assertIn("bad_clustal_format.aln", get_files)
        self.assertIn("no_ref_seq.aln", get_files)
        
    def test_alignio_hates_infile(self):
        """
            Since you're checking for a return code from sys.exit(), 
            make sure it's not confused with an exit code of 1 coming from 
                1) some other complaint raised by subprocess 
                2) some screwup in the tested script
        """

        checkit=subprocess.run(["python", "extract_demo.py", "-i", "bad_clustal_format.aln"], capture_output=True)
        #checkit=subprocess.run(["python", "../resource_files/run_hw.py", "-i", "../resource_files/bad_clustal_format.aln"], shell=True, text=True, capture_output=True)

        #print("\n\targs={}\n\treturn code={}\n\tstdout={}\n\tstderr={}".format(checkit.args, checkit.returncode, checkit.stdout, checkit.stderr))
        
        # in the future, if you need to search the error message 
        # for something like 'ValueError:' you can use regex and assertTrue()
        # but make sure that the stderr is a string by setting text=True
        #value_error=re.compile('ValueError:')
        #self.assertTrue(value_error.search(checkit.stderr))
        
        # try...except will suppress the error output.
        # good thing I used a return code
        self.assertEqual(checkit.returncode, 1)
        
    def test_no_ref_seq(self):
        """
            return code should be 1 if there's no E. coli reference sequence. 
            Either CP016007.2543965.2545520 or CP016007.3589827.3591382
        """
        
        checkit=subprocess.run(["python", "extract_demo.py", "-i", "no_ref_seq.aln"], capture_output=True)
        
        self.assertEqual(checkit.returncode, 1)
        
    def test_check_output(self):
        """
            the output goes to /testing/processed_data, which is nice
            
            probably not the best way to test, but since the extracted sequences should only be Ts, 
            then asserting that the set("TTTTT") is equal to the set("T"). Both should be just T.
        """
        
        checkit=subprocess.run(["python", "extract_demo.py", "-i", "good_clustal_format.aln", "-e", "5,13"], capture_output=True, text=True)
        
        print("\n\targs={}\n\treturn code={}\n\tstdout={}\n\tstderr={}".format(checkit.args, checkit.returncode, checkit.stdout, checkit.stderr))
        
        # processed_data/extracted_16s_regions_2019-10-27/extracted_16s_2019-10-27_2037.fna
        
        spl_output=checkit.stdout.strip().split("/")
        new_file="processed_data/{0}/{1}".format(spl_output[-2], spl_output[-1])
        
        for x in SeqIO.parse(new_file, "fasta"):
            self.assertEqual(set(x), set("T"))
          
        # clean up by removing the directory and all files in it
        shutil.rmtree("processed_data/{}".format(spl_output[-2]))
        
        


if __name__ == '__main__':
    unittest.main()