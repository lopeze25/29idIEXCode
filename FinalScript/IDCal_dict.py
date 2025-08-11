import ast
import datetime
import os
from os.path import join

import numpy as np
from datetime import date


"""

Writer, reader, and parser 
Some elements can be adjusted for readability, user experience, or ease of use—such as error handling and variable names.

The writer takes any entry and appends it to the specified file path.
The reader reads entries from a file. A major change is that it now uses tuples,
which improves indexing and searching by comment compared to my original idea of dictionaries which has unique keys
It can also handle flux curves if the format for Flux_Curves.txt is changed, and it can convert dictionaries of lists into proper dictionaries.

The parser finds the coefficient corresponding to a given entry.
Uses dictionaries
written by Ezekiel Lopez 20250627
"""

# testCal class to test similar to before 
global default_path, default_filename
#default_path=r"C:\Users\29iduser\Documents\GitHub\"
default_path = r"/home/beams/29IDUSER/Documents/User_Macros/Macros_29id/IEX_Dictionaries/"
#default_fname = "Dict_IDCal.txt"
default_filename = "Dict_IdCal_test.txt"

class IDCal_dict:
    def __init__(self):
        """
        written by Ezekiel Lopez 20250627
        Used to read and write ID calibration to a text file
        default filepath and filename are the global variables default_filename, default_path

        usage:
        reading in a dictionary and getting coefs
            cal = IDCal_dict()
            cal.read()
            cal.get_coefs
        updating the cal.new_entry 
            cal.update_entry
        writing to file the cal.new_entry
            cal.update_entry
        """
        self.last_comment = None
        self.last_entry = None
        self.new_comment = None
        self.new_entry = None

    def write_old(self, new_entry, **kwargs):
        """
        Appends/Writes a new calibration to the calibration file.
        Args: 
            new_entry (dict)
        kwargs: 
            filename (str)
            path (str)
            comment (str)
            debug (bool)
        """
        kwargs.setdefault('filename', default_filename)
        kwargs.setdefault('path', default_path)
        kwargs.setdefault('debug', False)
        kwargs.setdefault("new_entry_string","=======")
        kwargs.setdefault("comment",'')
    
        filename = kwargs['filename']
        debug = kwargs['debug']
    
        fpath = join(kwargs["path"],kwargs["filename"])
    
        # Make the directory if it doesnt exist
        os.makedirs(fpath, exist_ok=True)
    
        # Appends it to the file
        with open(fpath, "a+") as f:
    
            f.write("\n"+kwargs['new_entry_string']+" " + str(date.today()) + ": " + kwargs['comment'] + "\n")
            f.write(str(new_entry))
            f.write("\n")
    
        if debug:
            print('write_calibration fpath:', fpath)
            print('ID calibration:', new_entry)
        print("Success")
       
    def read(self, **kwargs):
        """
        Reads any dictionary entry from a file.
        writes self.last_entry = comment,entry
    
        Args:
            fpath (str): Path to the file to read from.
            kwargs:
                index (int): -1 to get the last entry, or specific index.
                comment (str): Search by comment string in comment.
                date (str): Search by date string in comment.
                debug (bool): If True, print debug information.
    
        """
        kwargs.setdefault('filename', default_filename)
        kwargs.setdefault('path', default_path)
        kwargs.setdefault("index", -1)
        kwargs.setdefault("comment", '')
        kwargs.setdefault("date", '')
        kwargs.setdefault("new_entry_string","=======")
        kwargs.setdefault("debug", False)
    
        fpath = join(kwargs["path"],kwargs["filename"])
    
        #create a list of all the entries to be able to find any enry init 
        entries = []
    
        with open(fpath, 'r') as f:
            content = f.read()
    
        #find the five ======, and it split 
    
        raw_entries = content.split(kwargs['new_entry_string']) #we had gotten rid of the comment indicator
        for e,raw in enumerate(raw_entries):
            #removes all leading and trailing whitespace
        
            if raw.strip():  
                    try:
                        #Find first { for the block
                        brace_index = raw.index("{")
                        comment = raw[:brace_index].strip()
                        dict_text = raw[brace_index:]
                        last_brace_index = dict_text.rfind("}")
                        clean_block = dict_text[:last_brace_index + 1]
                        
                        # Safely evaluate the dictionary
                        entry_dict = ast.literal_eval(clean_block)
                        
                        #create a tuple and add it to entries
                        entries.append((comment, clean_block))
                    except Exception as e:
                        if kwargs['debug']:
                            print("Error")
                            continue
        
        # Search entries from top to bottom for a match by comment or date
        index2 = 0
    
        while index2 < len(entries):
            comment, block = entries[index2]
            if (kwargs['comment']!='' and kwargs['comment'] in comment) or (kwargs["date"]!='' and kwargs["date"] in comment):
                kwargs["index"] = index2
                break
            index2 += 1
    
        index = kwargs['index']
        comment, entry = entries[index]
        
        """
        #Convert list-of-pairs to dict if needed
        for grating in entry:
            for mode in entry[grating]:
                if isinstance(entry[grating][mode], list):
                    entry[grating][mode] = { bkpt: coefs for bkpt, coefs in entry[grating][mode]}
    
        """
    
        if kwargs['debug']:
            print("\ Selected Entry comment", comment)
            print("\n Dictionary:", entry )
    
        self.last_comment = comment
        self.last_entry = ast.literal_eval(entry)
        # add p
             
    def get_coefs(self,entry, grt, IDmode, energy_eV):
        """
         parse entry only works for dictionaries of the type:
        {"MEG": {0: {bkpt1: [polynomials], bkpt2: [polynomials]}}}
    
        Parses a calibration dictionary and returns the coefficients corresponding
        to the first breakpoint greater than the specified energy value.
    
        Args:
            entry (dict): Calibration dictionary.
            grating (str): e.g., "MEG"
            ID_mode (int): e.g., 0
            energy_eV (float): photon energy in eV
    
        Returns:
            list: Coefficient list for the nearest breakpoint >= energy_eV
            bkpt_max = 3800
    
        """
        #check for grating
        if grt in entry.keys():
            if IDmode in entry[grt].keys():
                breakpoint_list = list(entry[grt][IDmode].keys())
                sorted_breakpoints = np.sort(np.array(breakpoint_list))
                #find breakpoint to use
                for bp in sorted_breakpoints:
                    if energy_eV < bp:
                        bkpnt = bp
                        coefs = entry[grt][IDmode][bkpnt]
                        
                    else:
                    #just returns the highest breakpoint if the eV is above all breakpoints 
                        bkpnt = sorted_breakpoints[-1]
                        coefs = entry[grt][IDmode][bkpnt]
            else:
                e = "IDmode not in entry"   
                coefs = [1]
        else:
            e = "grt not in entry"
            coefs = [1]
        return coefs
        
    def update_entry(self, grt, IDmode, bkpt, coefs,):
        """
        format for entry
        """
        d = {grt:{IDmode:{bkpt:coefs}}}
        if self.new_entry == None:
            self.new_entry = d
        else:
            grt = list(d.keys())[0]
            mode = list(d[grt].keys())[0]
            bkpt = list(d[grt][mode].keys())[0]
            
            if grt in self.new_entry.keys():
                if mode in self.new_entry[grt].keys():
                    #updating bkpnts
                    self.new_entry[grt][mode].update(d[grt][mode])
                else:
                    #update mode
                    self.new_entry[grt].update(d[grt])
            else:
                #update grt
                self.new_entry.update(d)          
    
    def write_new_entry(self,**kwargs):
        """
        writes to the dictionary the values in self.entry and self.comment 

        **kwargs:
            filename: default is global default_filename
            path: default is global default_path
            debug
        """
        kwargs.setdefault("new_entry_string","=======")
        
        kwargs.setdefault('filename', default_filename)
        kwargs.setdefault('path', default_path)
        kwargs.setdefault('debug', False)
        fpath = join(kwargs["path"],kwargs["filename"])
        
        #getting attributes to write
        if 'comment' in kwargs:
            self.comment = kwargs['comment']
        else:
            self.comment = (kwargs['new_entry_string'] + " " +str(date.today()) + ": " + kwargs['comment'])
            
        if self.new_entry == None:
            print('No entry defined; nothin written')
        
        #writing to dictionary file 
        else:         
            #os.makedirs(fpath, exist_ok=True)
            with open(fpath, "a+") as f:
                f.write("\n")
                f.write(str(self.new_comment) + "\n" + str(self.new_entry))
                f.write("\n")
            if kwargs['debug']: 
                print('write_calibration fpath:', fpath)
                print('ID calibration:', self.new_entry)