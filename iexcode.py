import ast
import datetime

import os

def write_calibration(new_entry, **kwargs):
    """
Appends/Writes a new calibration to the calibration file
Args: takes any dict(entry)
kwargs: 
filename (str)
path(str)
debug(bool) 
"""
    kwargs.setdefault('filename', "Dict_IDCal.txt")
    kwargs.setdefault('path', ".")
    kwargs.setdefault('debug', False)
    filename = kwargs.get('filename', "Dict_IDCal.txt")
    path = kwargs.get('path', ".")
    debug = kwargs.get('debug', False)

    fpath = join(path, filename)
    #Make sure directory is created
    os.makedirs(path, exist_ok=True)

    with open(fpath, "a+") as f:
        f.write("\n======= " + date.today() + ": " + kwargs.get("comment", "") + "\n")
        f.write(str(new_entry))
        f.write("\n")
    if kwargs['debug']:
        print('write_calibration fpath:',fpath)
        print('ID calibration:',new_entry)


filename = kwargs.get('filename', )
def write_calibration(new_entry, **kwargs):
    """
    Appends a new calibration entry to the calibration file.

    Args:
        new_entry (dict): Calibration dictionary.
        kwargs:
            filename (str): Output file name (default: IDcal_fname).
            path (str): Directory path (default: dictionary_path).
            debug (bool): If True, prints debug info.
    """
    filename = kwargs.get('filename', "Dict_IDCal.txt")
    path = kwargs.get('path', ".")
    debug = kwargs.get('debug', False)

    fpath = join(path, filename)
    os.makedirs(path, exist_ok=True)

    with open(fpath, "a+") as f:
        f.write("\n======= " + today() + ": " + kwargs.get("comment", "") + "\n")
        f.write(str(new_entry))
        f.write("\n")

    if debug:
        print(f"[DEBUG] Wrote to {fpath}")
        print("[DEBUG] Entry:\n", str(new_entry))



def read(self, fpath, **kwargs):
    """
    Reads any dictionary entry from a file .
    Args:
     fpath(str), Path to the desired file to read from
    kwargs:
    index, -1 to get the last entry
    comment,search by comment
    date,search by date
    debug
    Return desired entry
    """
    #Set the last entry by default
    kwargs.setdefault("index", -1)
    kwargs.setdefault("comment", None)
    kwargs.setdefault("date", None)
    kwargs.setdefault("debug", False)
    index = kwargs.get("index")
    comment = kwargs.get("comment")
    date = kwargs.get("date")
    debug = kwargs.get("debug")
    #List of all entries
    entries = []
    header = None
    block = []
  
    with open(fpath, 'r') as f:
        for line in f:
            if line.startswith("======="):
                if header is not None:
                    # Combine and strip lines from the block
                    combined_lines = ''
                    for l in block:
                        combined_lines += l.strip()
                    entries.append((header, combined_lines))
                #Stores the entire line into the header and remove whitespace 
                header = line.strip()
                block = []
            else:
                block.append(line)
# Make sure to add the last block after the loop
    if header and block:
        combined_lines = ''
        for l in block:
            combined_lines += l.strip()
        entries.append((header, combined_lines))

    # Search by comment or date in the header 
    i = len(entries) - 1
    while i >= 0:
        header, block = entries[i]
        if (comment in header) or (date in header):
            break
    i -= 1
    # Since we have all the entries in a list, find the index position of that list and reutnr hedad.blcok
    if index is not None:
        header, block = entries[index]
    #find the last_brace to clean up the block so it can evalvulate the string to dictionary 
    last_brace_index = block.rfind("}")
    clean_block = block[0:last_brace_index + 1]
    #Evalulate the string to a dictionary
    entry = ast.literal_eval(clean_block)

    return entry




import numpy as np

def get_coeff_list(entry, grating, ID_mode):
    return entry[grating][ID_mode]

def parse_entry(entry, grating, ID_mode, energy_eV):
  
    coeff_list = get_coeff_list(entry, grating, ID_mode)

    # Convert to numpy arrays for sorting
    breakpoints = np.array([pair[0] for pair in coeff_list])
    coefficients = np.array([pair[1] for pair in coeff_list], dtype=object)

    # Sort using numpy argsort
    sort_idx = np.argsort(breakpoints)
    breakpoints_sorted = breakpoints[sort_idx]
    coeffs_sorted = coefficients[sort_idx]

    # Find the first breakpoint greater than energy_eV
    for i, bp in enumerate(breakpoints_sorted):
        if energy_eV < bp:
            return coeffs_sorted[i].tolist()

    # If none found, return the last
    return coeffs_sorted[-1].tolist()


def parse_entry(entry, grating, ID_mode, energy_eV):
    """
    Parses a calibration dictionary and returns the coefficients corresponding
    to the first breakpoint greater than the specified energy value.

    Args:
        entry (dict): Calibration dictionary of the form {grating: {ID_mode: {breakpoint: coeffs}}}
        grating (str): e.g., "MEG"
        ID_mode (str): e.g., "RCP"
        energy_eV (float): photon energy in eV

    Returns:
        list: Coefficient list for the nearest breakpoint ≥ energy_eV
    """

    # Extract the relevant sub-dictionary for the given grating and mode
    breakpointdict = entry[grating][ID_mode]

    # Extract breakpoints from the dictionary keys (integers or floats)
    breakpoint_list = list(breakpointdict.keys())  # This gives a flat list

    # Sort the breakpoints using NumPy
    sorted_breakpoints = np.sort(np.array(breakpoint_list))

    # Loop over the sorted list to find the first breakpoint larger than energy_eV
    for bp in sorted_breakpoints:
        if energy_eV < bp:
            return breakpointdict[bp]

    # If energy is greater than all breakpoints, return the highest one
    return breakpointdict[sorted_breakpoints[-1]]
