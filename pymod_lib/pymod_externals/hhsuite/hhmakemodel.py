#!/usr/bin/env python

from pymod_lib.pymod_externals.hhsuite.hh_reader import read_result
from copy import deepcopy
import re, os, sys, tempfile, glob

from operator import itemgetter  # hzhu
from itertools import groupby    # hzhu

EMPTY = '*'
GAP = '-'
DEBUG_MODE = False

class Gap: 
    """ A gap is a continuous stretch of indels.
        It is defined by a opening position and a size/length
    """
    def __init__(self, open_pos, size):
        self.open_pos = open_pos   # gap opening position
        self.size = size           # num of indels in the gap

    def __repr__(self):
        return 'Gap opening pos = %d, size = %d' % (self.open_pos, self.size)
        
class Grid:
    """
    Implementation of 2D grid of cells
    Includes boundary handling
    """
    
    def __init__(self, grid_height, grid_width):
        """
        Initializes grid to be empty, take height and width of grid as parameters
        Indexed by rows (left to right), then by columns (top to bottom)
        """
        
        self._grid_height = grid_height
        self._grid_width = grid_width
        self._cells = [ [ EMPTY for dummy_col in range(self._grid_width) ]
                       for dummy_row in range(self._grid_height)]
                
    def __str__(self):
        """ Return multi-line string represenation for grid """
        
        ans = ''
        for row in range(self._grid_height):
            ans += ''.join(self._cells[row])
            ans += '\n'
        return ans

    def clear(self):
        """ Clears grid to be empty """
        
        self._cells = [[EMPTY for dummy_col in range(self._grid_width)]
                       for dummy_row in range(self._grid_height)]
    
    def get_grid_height(self):
        """ Return the height of the grid """
        
        return self._grid_height

    def get_grid_width(self):
        """ Return the width of the grid """
        
        return self._grid_width
    
    def get_cell(self, row, col):
        return self._cells[row][col]

    def get_seq_start(self, row):
        """ Returns the start position of the sequence """
        
        index = 0
        for pos in self._cells[row]:
            if pos != EMPTY: 
                return index
            index += 1

        return None

    def get_seq_end(self, row):
        """ Returns the end position of the sequence """
        
        index = 0
        for pos in reversed(self._cells[row]):
            if pos != EMPTY: 
                return self.get_grid_width() - index
            index += 1

        return None

    def get_gaps(self, row): 
        """ Return the position of gaps in a row """
        
        gaps = list()

        index = 0
        for pos in self._cells[row]:
            if pos == GAP: 
                gaps.append(index)
            index += 1

        return gaps

    def get_gaps_ref_gapless(self, row):
        """ Return the pos of gaps in a row.
            The opening positions of the gaps are wrt. the gapless seq
        """
        # get all the indels
        indels = self.get_gaps(row)  
        gaps = []
        # combine continuous indels into a gap
        for k,i in groupby( enumerate(indels), lambda x: x[0]-x[1] ):
            g = list(map(itemgetter(1), i))
            gaps.append( Gap(g[0], len(g)) )

        # offset the gap opening positions
        for i in range(1, len(gaps)):
            # offset by total gap number before
            gaps[i].open_pos -= sum([gaps[j].size for j in range(i)])
            
        return gaps  # a list of Gap instances
    
    def get_seq_indeces(self, row):

        seq = list()
        for pos, res in enumerate(self._cells[row]):
            if res != EMPTY and res != GAP:
                seq.append(pos)

        return seq
            
    ## def get_gap_list(self):  # hzhu commented this out. wrote a new version
    ##     """ Returns a list of list of all gap positions in the sequence grid. """
    ##     gap_pos = set()

    ##     for row in range(self.get_grid_height()):
    ##         for gap in self.get_gaps(row):
    ##             gap_pos.add(gap)

    ##     gap_pos = list(sorted(gap_pos))
    ##     boundaries = [ (x + 1) for x, y in zip(gap_pos, gap_pos[1:]) if y - x != 1 ]
 
    ##     gap_list = list()
    ##     prev = 0
 
    ##     for boundary in boundaries:
    ##         sub_list = [ pos for pos in gap_pos[prev:] if pos < boundary ]
    ##         gap_list.append(sub_list)
    ##         prev += len(sub_list)
 
    ##     gap_list.append([ x for x in gap_pos[prev:]])
    
    ##     return gap_list

    def get_gap_list(self):
        """ Returns a list of Gap instances for all rows in the grid
        """
        gap_dict = dict()  # each position should occur as gap at most once
                           # keys are gap openning positions
                           # values are Gap instances
        gap_list = []
        for row in range(self.get_grid_height()):
            gap_pos = []
            gaps = self.get_gaps_ref_gapless(row)
             
            for g in gaps:
                if g.open_pos in gap_dict: # if there is already gaps at this open pos
                    if g.size > gap_dict[g.open_pos].size: # if new gap is bigger
                        gap_dict[g.open_pos] = g  # keep the larger gap as they overlap
                else:
                    gap_dict[g.open_pos] = g
                    
        gap_list = sorted(list(gap_dict.values()), key=lambda x: x.open_pos) # sort according to start position
        return gap_list  # a list of Gap instances
        
    def set_gap(self, row, col):
        """ Set cell with index (row, col) to be a gap """
        
        self._cells[row][col] = GAP

    def set_empty(self, row, col):
        """ Set cell with index (row, col) to be a gap """
        
        self._cells[row][col] = EMPTY
    
    def set_cell(self, row, col, res):
        """ Set cell with index (row, col) to be full """
        
        self._cells[row][col] = res

    def is_empty(self, row, col):
        """ Checks whether cell with index (row, col) is empty """
        
        # return self._cells[row][col] == EMPTY
        try:
            return self._cells[row][col] == EMPTY
        except IndexError:
            print("WARNING!")
            return True
        

    def is_gap(self, row, col):
        """ Checks whetehr cell with indxex (row, col) is a gap """
        
        return self._cells[row][col] == GAP

    def insert_gaps(self, cols):
        """ Inserts a gaps into a column of the template grid """
        
        for col in cols:
            for row in range(self._grid_height):
                if col >= self.get_seq_start(row) and col < self.get_seq_end(row):
                    self._cells[row].insert(col, GAP)
                else:
                    self._cells[row].insert(col, EMPTY)
            
            self._grid_width += 1

    def insert_gaps_row(self, cols, row):
        """ Intert gaps into cols only for certain row"""
        for col in cols:
            if col >= self.get_seq_start(row) and col < self.get_seq_end(row):
                self._cells[row].insert(col, GAP)
            else:
                self._cells[row].insert(col, EMPTY)
            # NOTE: grid_with should not be changed after every row is updated.    
            #self._grid_width += 1

    def clean_trail_empty(self):
        """ Remove all trailing EMPTY and pad grid to same width"""
        # first find out the max length (exluding trailing EMPTY)
        max_width = 0
        for row in range(self._grid_height):
            for i in range(len(self._cells[row])-1, -1, -1):
                if self._cells[row][i] != EMPTY:
                    break
            if i+1 > max_width:
                max_width = i+1
                
        # delete excessive EMPTY        
        for row in range(self._grid_height):
            del self._cells[row][max_width:]

        # then pad all rows to the same length
        [self._cells[row].append( EMPTY * (max_width-len(self._cells[row])) ) \
                                 for row in range(self._grid_height) if len(self._cells[row]) < max_width]
        self._grid_width = max_width
        return
    
    def remove_gaps(self, keep_width=True): # hzhu add keep_width option
        """ Removes all gaps from the grid. """

        for row in range(self.get_grid_height()):
            not_gap = list()
            for col in range(self.get_grid_width()):
                if not self.is_gap(row, col):
                    not_gap.append(col)

            self._cells[row] = [ self._cells[row][col] for col in not_gap ]

            if keep_width:  # hzhu only pad to original width if desired
                for del_pos in range(self._grid_width - len(not_gap)):
                    self._cells[row].append(EMPTY)

        if not keep_width: # hzhu if width is not kept, make sure width is consistent
            self.clean_trail_empty()
                    
        return
                

class QueryGrid(Grid):

    def __init__(self, grid_height, grid_width):
        Grid.__init__(self, grid_height, grid_width)
    
    def get_query_start(self, row):
        """ Returns the query start position """
        return self.get_seq_start(row) + 1

    def get_query_end(self, row):
        """ Returns the query end postion """
        
        return self.get_seq_end(row) - len(self.get_gaps(row))

    def get_col_residue(self, col):
        """ Tries to find a the query residue in a given column. Used by derive_global_seq() to
        identify the global query sequence """
        
        for row in range(self.get_grid_height()):
            if not self.is_empty(row, col):
                return self._cells[row][col]

        return GAP

class TemplateGrid(Grid):

    def __init__(self, grid_height, grid_width):
        Grid.__init__(self, grid_height, grid_width)

        self._start = list()
        self._end = list()
        self._pdb_code = list()
        self._chain = list()
        self._organism = list()
        self._resolution = list()
        
    def display(self):
        """ Return multi-line string represenation for grid """
        
        ans = ''
        for row in range(self._grid_height):
            ans += '>P1;{p}\nstructure:{p}:{s}:{c}:{e}:{c}::{o}:{r}:\n{a}*\n'.format(
                p = self._pdb_code[row],
                s = add_white_space_end(self.get_template_start(row), 4),
                e = add_white_space_end(self.get_template_end(row), 4),
                c = self._chain[row],
                o = self._organism[row],
                r = self._resolution[row], 
                a = ''.join(self._cells[row]).replace(EMPTY, GAP).replace('#', GAP))

        return ans

    def debug(self, row):
        """ Return multi-line string represenation for grid, for debugging purposes """

        ans = '{p}\nInternal: {s}, {e} Query: {qs}, {qe} Gaps ({g1}): {g2}\n{seq}\n'.format(
            p = self._pdb_code[row],
            s = self.get_seq_start(row),
            e = self.get_seq_end(row),
            qs = self.get_template_start(row),
            qe = self.get_template_end(row),
            g1 = len(self.get_gaps(row)),
            g2 = ', '.join([str(gap) for gap in self.get_gaps(row)]),
            seq = ''.join(self._cells[row]))

        return ans 

    def set_metadata(self, row, start, end, pdb_code, chain, organism, resolution):
        """ Used by create_template_grid() to setup metadata of pir template """
        
        self._start.append(start)
        self._end.append(end)
        self._pdb_code.append(pdb_code)
        self._chain.append(chain)
        self._organism.append(organism)
        self._resolution.append(resolution)

    def set_map(self, row, start, end):
        
        self._start[row] = start
        self._end[row] = end

    def get_template_start(self, row):
        """ Returns the template start position """
        
        return self._start[row]

    def get_template_end(self, row):
        """ Return sthe template end position """
        
        return self._end[row]

    def del_row(self, row):
        """ Removes a complete template entry from the grid """
        
        del self._cells[row]
        del self._start[row]
        del self._end[row]
        del self._pdb_code[row]
        del self._chain[row]
        del self._organism[row]
        del self._resolution[row]
        self._grid_height -= 1

# Helper functions

def add_white_space_end(string, length):
    """ Adds whitespaces to a string until it has the wished length"""
    
    edited_string = str(string)
    
    if len(edited_string) >= length:
        return string
    else:
        while len(edited_string) != length:
            edited_string += ' '
    
    return edited_string


def get_query_name(hhr_file):

    with open(hhr_file) as fh:
        for line in fh:
            if line.startswith('Query'):
                # match the PDB Code
                m = re.search('(\d[A-Z0-9]{3})_(\S)', line) 

                if m: 
                    pdb_code = m.group(1)
                    chain = m.group(2)
                else: 
                    pdb_code = 'UKNP'
                    chain = 'A'
                    # raise ValueError('Input HHR-File Does not seem to be a PDB-Structure')

                break

    return pdb_code, chain


def template_id_to_pdb(template_id):
    """
    Extracts PDB ID and chain name from the provided template id
    """
    # match PDBID without chain (8fab, 1a01)
    m = re.match(r'/^(\d[A-Za-z0-9]{3})$', template_id)
    if m:
        return m.group(1).upper(), 'A'
    
    # PDB CODE with chain Identifier
    m = re.match(r'^(\d[A-Za-z0-9]{3})_(\S)$', template_id)
    if m:
        return m.group(1).upper(), m.group(2).upper()

    # Match DALI ID
    m = re.match(r'^(\d[A-Za-z0-9]{3})([A-Za-z0-9]?)_\d+$', template_id)
    if m:
        return m.group(1).upper(), m.group(2).upper()
    
    # No PDB code and chain identified
    return None, None


def create_template_grid(hhr_data):
    """ Creates a template grid """

    total_seq = len(hhr_data)
    templ_max = max( [ hhr.start[0] + len(to_seq(hhr.template_ali)) for hhr in hhr_data ] ) - 1


    template_grid = TemplateGrid(total_seq, templ_max)

    for row, template in enumerate(hhr_data):
        seq_start = template.start[0] - 1
        templatealignment = to_seq(template.template_ali)
        seq_end = seq_start + len(templatealignment)

        # Load Meta Data
        start = template.start[1]
        end = template.end[1]

        # Get pdb_code and chain identifier of template
        pdb_code, chain =  template_id_to_pdb(template.template_id)

        m = re.search("(\d+.\d+)A", template.template_info) # try to extract resolution of the structure
        
        if m: 
            resolution = m.group(1)
        else: 
            resolution = ""

        m = re.search("\{(.*)\}", template.template_info) # try to extract the organism
        if m: 
            organism = m.group(1).replace(":", " ") # make sure that no colons are in the organism
        else: 
            organism = ""

        template_grid.set_metadata(row, start, end, pdb_code, chain, organism, resolution)

        # Write sequence into the grid
        for pos, col in enumerate(range(seq_start, seq_end)):
            template_grid.set_cell(row, col, templatealignment[pos])

    return template_grid


def to_seq(ali):
    if isinstance(ali, list):
        return ''.join(ali)
    else:
        return ali
    

def create_query_grid(hhr_data):
    """ Creates a Query Grid """
    
    total_seq = len(hhr_data)
    query_max = max( [ hhr.start[0] + len(to_seq(hhr.query_ali)) for hhr in hhr_data ] ) - 1

    query_grid = QueryGrid(total_seq, query_max)

    for row, query in enumerate(hhr_data):

        queryalignment = to_seq(query.query_ali)
        query_start = query.start[0] - 1
        query_end = query_start + len(queryalignment)

        for pos, col in enumerate(range(query_start, query_end)):
            if queryalignment[pos] not in ['Z', 'U', 'O', 'J', 'X', 'B']: # CAUTION

                query_grid.set_cell(row, col, queryalignment[pos])

    return query_grid

def create_gapless_grid(grid):
    """ Returns a gapless grid """

    gapless = deepcopy(grid)
    gapless.remove_gaps(keep_width=False)  # hzhu:  shrink grid

    return gapless

def process_query_grid(query_grid, gapless_grid):
    """ Processes a query grid sucht that it contains all gaps
    """
    gaplist = query_grid.get_gap_list()
    off_set = 0
    
    for g in gaplist:
        gapless_grid.insert_gaps([ p + off_set for p in range(g.open_pos, g.open_pos+g.size) ])
        off_set += g.size
        
    return gapless_grid

def derive_global_seq(processed_query_grid, query_name, query_chain):

    global_seq = list()

    for col in range(processed_query_grid.get_grid_width()):
        global_seq.append(processed_query_grid.get_col_residue(col))

    # this is the query entry
    header = '>P1;{q}\nsequence:{q}:1    :{c}:{l}  :{c}::::\n'.format(
        q = query_name, 
        l = len(global_seq),
        c = query_chain)

    return header + ''.join(global_seq) + '*'

def process_template_grid(query_grid, template_grid):
    """ Insertes Gaps into the template grid
        Only add gaps from **other** query_grids into template grid (NOT gapless)
    """
    gaplist = query_grid.get_gap_list()  # use this to keep the offset
    
    for row in range(template_grid.get_grid_height()):
        # do NOT consider gaps in current query row         
        gaplist_row = query_grid.get_gaps_ref_gapless(row)
        gapdict_row = dict(zip([g.open_pos for g in gaplist_row],
                               [g.size     for g in gaplist_row]))
        off_set = 0
        for g in gaplist:
            # if there is a gap with same opening position in the current row,
            # only consider g if it is larger than the on in the current row
            if g.open_pos in gapdict_row:
                if g.size > gapdict_row[g.open_pos]: 
                    template_grid.insert_gaps_row([ p + off_set for p in range(g.open_pos,
                                                                              g.open_pos+g.size-gapdict_row[g.open_pos]) ], row)
            else:
                template_grid.insert_gaps_row([ p + off_set for p in range(g.open_pos, g.open_pos+g.size) ], row)
                 
            off_set += g.size  # even if the gaps are not inserted, the offset should be adjusted

    template_grid.clean_trail_empty()  # clean the redundant trailing EMPTY char
    
    return template_grid


def remove_self_alignment(template_grid, query_name):
    """ Removes a self alignment from the final pir alignment to prevent clashes with MODELLER """
    
    to_delete = list()
    
    for row in range(template_grid.get_grid_height()):
        if template_grid._pdb_code[row] == query_name:
            to_delete.append(row)

    for row in reversed(to_delete):
        template_grid.del_row(row)

    return True

def write_to_file(line_list, fname):
    """ Writes the final pir file """
    
    with open(fname, 'w+') as fout:
        for line in line_list:
            fout.write(line + "\n")

def arg():
    import argparse
    description = """Creates a MODELLER alignment (*.pir) from a HHSearch results file (*.hhr)."""
    epilog= '2016 Harald Voehringer.'
    # Initiate a ArgumentParser Class
    parser = argparse.ArgumentParser(description = description, epilog = epilog)
    
    # Call add_options to the parser
    parser.add_argument('input', help = 'results file from HHsearch with hit list and alignment', metavar = 'FILE')    
    parser.add_argument('cifs', help = 'path to the folder containing cif files', metavar = 'DIR')
    parser.add_argument('pir', help = 'output file (PIR-formatted multiple alignment)', metavar = 'FILE')
    parser.add_argument('output', help = 'path to the folder where modified cif files should be written to', metavar = 'DIR')

    parser.add_argument('-v', '--verbose', action = 'store_true', help = 'verbose mode')
    parser.add_argument('-m', nargs = '+', help = 'pick hits with specified indices (e.g. -m 2 5)', metavar = 'INT')
    parser.add_argument('-e', type = float, help = 'maximum E-Value threshold (e.g. -e 0.001)', metavar = 'FLOAT')
    parser.add_argument('-r', type = float, help = 'residue ratio (filter alignments that have contribute at least residues according to the specified ratio).',
        default = 0, metavar = 'FLOAT')

    parser.add_argument('-c', help = 'convert non-canonical residues (default = True)', action = 'store_true', default = True)
    

    return parser


class PyMod_args:
    def __init__(self, args_dict):
        for arg in args_dict:
            setattr(self, arg, args_dict[arg])
        print(self.input)


def main(pymod_args={}):
    
    if not pymod_args:
        import sys
        parser = arg()
        args = parser.parse_args(sys.argv[1:])
    else:
        args = PyMod_args(pymod_args)
    
    global DEBUG_MODE
    
    if args.verbose:
        DEBUG_MODE = True

    query_name, query_chain = get_query_name(args.input)

    data = read_result(args.input)
    selected_templates = list()

    if args.m and not args.e:
        selection = map(lambda x: int(x), args.m)
        print ('Selected templates {st}.'.format(st = ', '.join(args.m)))
        
        for i in selection:
            tmp_info = str(data[i - 1].template_info.split('>')[1])
            print ('{i}: {t}'.format(
                i = i,
                t = tmp_info[0:80]))

            selected_templates.append(data[i - 1])

    elif args.e and not args.m:
        print ('Selected templates satisfying E-val <= {e}'.format(e = args.e))
        
        e_values = { float(j.evalue):i for i, j in enumerate(data) }
        selection = sorted([ val for key, val in e_values.items() if key <= args.e ])

        for i in selection:
            tmp_info = str(data[i - 1].template_info.split('>')[1])
            print ('{i}: {t}'.format(
                i = i + 1,
                t = tmp_info[0:80]))

            selected_templates.append(data[i - 1])

    elif args.m and args.e:
        print ('! Please do not use option -m and -e at the same time ! Exiting.')
        sys.exit()
    else:
        selected_templates = data
        
        print ('Creating pir file using all templates ({n})'.format(
            n = len(selected_templates)))

    query_grid = create_query_grid(selected_templates) # load query grid
    print ('query_grid')
    print(query_grid)
    gapless_query_grid = create_gapless_grid(query_grid) # remove gaps
    print ('gapless_query_grid')
    print(gapless_query_grid)
    processed_query_grid = process_query_grid(query_grid, gapless_query_grid) # insert gaps
    ##processed_query_grid = process_query_grid(query_grid, query_grid) # insert gaps
    print ('processed_query_grid')
    print (processed_query_grid)
    glob_seq = derive_global_seq(processed_query_grid, query_name, query_chain) # derive query sequence
    template_grid = create_template_grid(selected_templates) # create template grid
    print ('template_grid')
    print (template_grid)
    processed_template_grid = process_template_grid(query_grid, template_grid) # insert gaps to template sequnces
    print ('processed_query_grid')
    print (processed_query_grid)
    print ('hzhu processed_template_grid')
    print (processed_template_grid)
    # final_grid = compare_with_cifs(processed_template_grid, args.cifs, args.output, args.c, args.r) # compare with atom section of cifs
    final_grid = processed_template_grid
    remove_self_alignment(final_grid, query_name) # remove self alignment if any
    write_to_file([glob_seq, final_grid.display()], args.pir)


if __name__ == "__main__":
    main()
