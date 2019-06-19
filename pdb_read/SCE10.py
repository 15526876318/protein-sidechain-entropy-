#!usr/bin/env python

from __future__ import print_function
import numpy
import re
import gzip
from Bio.PDB.PDBExceptions import PDBConstructionException
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.StructureBuilder import StructureBuilder

class PDBP_read(object):
	def __init__(self, get_header=False, structure_builder=None,PERMISSIVE=True):
		"""arguments:
			PERMISSIVE, Evaluated as a Boolean. If ture, the exception are caught, 
			some residues or atoms will be missing.THESE EXCEPTIONS ARE DUE TO PROBLEMS IN THE PDB FILE!
			structure_builder, an optional user implemented StructureBuilder class.
			
			"""
		#get a structure_builder class 
		if structure_builder is not None:
			self.structure_builder = structure_builder
		else:
			self.structure_builder = StructureBuilder()
		self.header = None
		self.trailer = None
		self.line_counter = 0
		self.PERMISSIVE = bool(PERMISSIVE)
		

	def get_structure(self,id,file):
		"""return the structure.
		   argurements:
			-id - the name of the sturecture
			-file - pdb filename
		"""
		self.header = None
		self.trailer = None
	#make a StructureBuilder instance
		self.structure_builder.init_structure(id)
		if file[-2:] == 'gz':
			#try:
				#with open(file,'r+',encoding='utf-8') as handle:
					#self._parse(handle.readlines())
			with gzip.open(file, 'r') as handle:   #######按照行读取pdb文件
				self._parse(handle.read().decode("utf-8").split('\n'))
		else:
			try:
				with open(file, 'r') as handle:
					self._parse(handle.readlines())
			except Exception:
				print("%s cannot be open!" % (file))
				#exit()
		#???????
		self.structure_builder.set_header(self.header)
		# return the structure instance
		structure = self.structure_builder.get_structure()
		return structure
		
	def get_header(self):
		"""return the header"""
		return self.header

	def get_trailer(self):
		"""return the trailer"""
		return self.trailer

	#private methods
	def _parse(self, header_coords_trailer):
		"""parser the pdb file(private)"""
		self.coords_trailer = self._get_header(header_coords_trailer)
		## parse the atomicdata; return the pdb file triler
		self.trailer = self._parse_coordinates(self.coords_trailer)

	def _get_header(self,header_coords_trailer):
		"""get the header of the pdb file"""
		structure_builder = self.structure_builder
		i = 0
		line_nums = len(header_coords_trailer)
		'''
		for index, line in enumerate(header_coords_trailer):
			print(index, line)
		print(line_nums, type(header_coords_trailer))
		'''
		for i in range(0,line_nums):
			structure_builder.set_line_counter(i + 1)
			line = header_coords_trailer[i]
			record_type = line[0:6]
			if record_type == "ATOM  " or record_type == "HETATM" or record_type == "MODEL ":
				break
		#header = header_coords_trailer[0:i]
		#return the rest of the coodstrailer
		self.line_counter = i
		coords_trailer = header_coords_trailer[i:]
		#header_dict = self._parse_pdb_header_list(header)
		return coords_trailer

	def _parse_coordinates(self,coords_trailer):
		"""parse the atomic data in teh PDB file """
		local_line_counter = 0
		#n=0
		structure_builder = self.structure_builder
		current_model_id = 0
		# Flag we have an open model
		model_open = 0
		current_chain_id = None
		current_segid = None
		current_residue_id = None
		current_resname = None
		lines_num1 = len(coords_trailer) 
		for i in range(0, lines_num1):
			line = coords_trailer[i].rstrip('\n')
			record_type = line[0:6]
			global_line_counter = self.line_counter + local_line_counter + 1
			# the all lines nums include header coods and trailer
			structure_builder.set_line_counter(global_line_counter)
			if record_type == "ATOM  ": #or record_type == "HETATM":
				#Initialize the Model - there was no explicit MODEL record
				if not model_open:
					structure_builder.init_model(current_model_id)
					current_model_id +=1
					model_open = 1
				fullname = line[12:16]
				# get rid of whitespace in atom names
				split_list = fullname.split()
				if len(split_list) !=1:
					# a atom has several species, eg "N B"
					name = fullname
				else:
					#eg: "CA"
					name = split_list[0]
				altloc = line[16]
				resname = line[17:20]
				chainid = line[21]
				try:
					serial_number = int(line[6:11])
				except Exception:
					serial_number = 0
				resseq = int(line[22:26].split()[0])
				icode = line[26]
				if record_type == "HETATM":
					if resname == "HOH" or resname == "WAT":
						hetero_flag = "W"
					else:
						hetero_flag = "H"
				else:
					hetero_flag = " "
				residue_id = (hetero_flag, resseq, icode)
				try:
					x = float(line[30:38])
					y = float(line[38:46])
					z = float(line[46:54])
				except Exception:
					raise PDBConstructionException("Invalid or missing coordinate(s) at line %i." % global_line_counter)
				coord = numpy.array((x, y, z), "f")
				try:
					occupancy = float(line[54:60])
				except Exception:
					self._handle_PDB_exception("Invalid or missing occupancy", global_line_counter)
					#occupancy = None  # Rather than arbitrary zero or one
					occupancy = None  # Rather than arbitrary zero or one
				if occupancy is not None and occupancy < 0:
					warnings.warn("Negative occupancy in one or more atoms", PDBConstructionWarning)
				try:
					bfactor = float(line[60:66])
				except Exception:
					self._handle_PDB_exception("Invalid or missing B factor",
											   global_line_counter)
					bfactor = 0.0  # The PDB use a default of zero if the data is missing
				segid  = line[72:76]
				element = line[76:78].strip().upper()
				if current_segid != segid:
					current_segid = segid
					structure_builder.init_seg(current_segid)
				if current_chain_id != chainid:
					current_chain_id = chainid
					structure_builder.init_chain(current_chain_id)
					current_residue_id = residue_id
					current_resname = resname
					try:
						structure_builder.init_residue(resname, hetero_flag, resseq, icode)
					except PDBConstructionException as message:
						self._handle_PDB_exception(message, global_line_counter)
				elif current_residue_id != residue_id or current_resname != resname:
					current_residue_id = residue_id
					current_resname = resname
					try:
						structure_builder.init_residue(resname, hetero_flag, resseq, icode)
					except PDBConstructionException as message:
						self._handle_PDB_exception(message, global_line_counter)
				# init atom
				try:
					structure_builder.init_atom(name, coord, bfactor, occupancy, altloc,
												fullname, serial_number, element)
				except PDBConstructionException as message:
					self._handle_PDB_exception(message, global_line_counter)
			elif record_type == "MODEL ":
				try:
					serial_num = int(line[10:14])
				except Exception:
					self._handle_PDB_exception("Invalid or missing model serial number",
											   global_line_counter)
					serial_num = 0
				structure_builder.init_model(current_model_id, serial_num)
				current_model_id += 1
				model_open = 1
				current_chain_id = None
				current_residue_id = None
			elif record_type == "END   " or record_type == "CONECT":
				# End of atomic data, return the trailer
				self.line_counter += local_line_counter
				return coords_trailer[local_line_counter:]
			elif record_type == "ENDMDL":
				model_open = 0
				current_chain_id = None
				current_residue_id = None
			elif record_type == "SIGUIJ":
				# standard deviation of anisotropic B factor
				siguij = [float(x) for x in (line[28:35], line[35:42], line[42:49],
											 line[49:56], line[56:63], line[63:70])]
				# U sigma's are scaled by 10^4
				siguij_array = (numpy.array(siguij, "f") / 10000.0).astype("f")
				structure_builder.set_siguij(siguij_array)
			elif record_type == "SIGATM":
				# standard deviation of atomic positions
				sigatm = [float(x) for x in (line[30:38], line[38:45], line[46:54],
											 line[54:60], line[60:66])]
				sigatm_array = numpy.array(sigatm, "f")
				structure_builder.set_sigatm(sigatm_array)
			local_line_counter += 1
		# EOF (does not end in END or CONECT)
		self.line_counter = self.line_counter + local_line_counter
		return []
				#info = (resname, resseq, serial_number,fullname, coord)
				#yield info

	def _handle_PDB_exception(self, message, line_counter):
		message = "%s at line %i." % (message, line_counter)
		if self.PERMISSIVE:
			# just print a warning - some residues/atoms may be missing
			warnings.warn("PDBConstructionException: %s\n"
						  "Exception ignored.\n"
						  "Some atoms or residues may be missing in the data structure."
						  % message, PDBConstructionWarning)
		else:
			# exceptions are fatal - raise again with new message (including line nr)
			raise PDBConstructionException(message)

if __name__ == "__main__":
	import sys
	p = PDBP_read(PERMISSIVE=True)
	filename = sys.argv[1]
	s = p.get_structure("scr", filename)




