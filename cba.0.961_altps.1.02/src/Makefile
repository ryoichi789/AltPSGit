CPP	= g++
#CPPFLAGS = -O2 -march=athlon64
CPPFLAGS = -O2 -march=core2
#CPPFLAGS = -O2
AR = ar
ARFLAGS = -rv
.SUFFIXES: .o .cpp .h
.cpp.o:
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $<

LIBDIR = ../lib
INCDIR = ../include
BINDIR = ../bin
LIBCBA = $(LIBDIR)/libcba.a
LIBMAIN = $(LIBDIR)/libcbamain.a
LIBALL = $(LIBMAIN) $(LIBCBA)

BOOST_LIBDIR = /usr/local/lib # for 4.1
BOOST_INCDIR = /usr/local/include/boost-1_31 # for 4.1

INCLUDE = # to use the Boost library, comment out this line and make the next line effective
#INCLUDE = -I$(BOOST_INCDIR) # for 4.1
LIBS = -lm -L$(LIBDIR) -lcbamain -lcba # to use the Boost library, comment out this line and make the next line effective
#LIBS = -lm -L$(LIBDIR) -lcbamain -lcba -L$(BOOST_LIBDIR) -lboost_regex-gcc-1_31 # for 4.1
LDFLAGS =

HEADERS = app.h exception.h\
	bio_sequence.h fasta_reader.h bio_seq_cmp.h\
	bio_seq_cmp_identical.h bio_seq_equiv_region.h\
	dna_sequence.h genetic_code.h\
	bio_seq_score_matrix.h bio_seq_alignment.h bio_seq_cmp_dp.h\
	bio_seq_db_search.h\
	region.h\
	hydropathy.h trans_membrane_prediction.h\
	pdb_record.h pdb_entry.h pdb_reader.h\
	atom.h position.h matrix.h molecule.h protein.h\
	sdf_entry.h sdf_reader.h\
	mol2_entry.h mol2_reader.h\
	secondary_structure.h asa_calculator.h\
	superposer.h structure_comparison.h\
	bit_string.h my_math.h cba_entry.h\
	dir_scan.h pdb_file_processor.h\
	ncbi_nrpdb.h bond_maker.h\
	date_time.h\
	mol_fin.h\
	use_boost.h\
	locarea.h altpsm.h

OBJS = \
	bio_sequence.o fasta_reader.o bio_seq_cmp.o bio_seq_cmp_identical.o\
	dna_sequence.o genetic_code.o\
	bio_seq_score_matrix.o bio_seq_alignment.o bio_seq_cmp_dp.o\
	bio_seq_db_search.o\
	hydropathy.o trans_membrane_prediction.o\
	pdb_record.o pdb_entry.o pdb_reader.o\
	atom.o position.o matrix.o molecule.o protein.o\
	sdf_entry.o sdf_reader.o\
	mol2_entry.o mol2_reader.o\
	secondary_structure.o asa_calculator.o\
	superposer.o structure_comparison.o\
	bit_string.o my_math.o cba_entry.o\
	dir_scan.o pdb_file_processor.o\
	ncbi_nrpdb.o bond_maker.o\
	mol_fin.o\
	date_time.o\
	locarea.o altpsm.o

APP_MAIN = app.o main.o

APPS = app_sample app_tmpl\
	compare_genomes find_orfs\
	align_seq search_seqdb\
	find_motif predict_tm\
	sum_pdb show_pdb_molinfo\
	sum_sdf show_atom_types\
	define_secondary_structures calc_asas\
	test_superposer compare_structures\
	pdb_to_cba count_seqs\
	local_time find_identical_mols morgan_numbering\
	get_pdb_atom_aas mol_fin_make mol_fin_cmp\
	altps

all: lib $(APPS)
lib: $(LIBALL)

######################################################################
# applications
######################################################################

# 1.2
app_sample: app_sample.o $(OBJS)
	$(CPP) -o $@ app_sample.o $(LIBS)
app_sample.o: app_sample.cpp app_sample.h $(HEADERS)

app_tmpl: app_tmpl.o $(OBJS)
	$(CPP) -o $@ app_tmpl.o $(LIBS)
app_tmpl.o: app_tmpl.cpp app_tmpl.h $(HEADERS)

# 2.1
compare_genomes: compare_genomes.o $(OBJS)
	$(CPP) -o $@ compare_genomes.o $(LIBS)
compare_genomes.o: compare_genomes.cpp compare_genomes.h $(HEADERS)

# 2.2
find_orfs: find_orfs.o $(OBJS)
	$(CPP) -o $@ find_orfs.o $(LIBS)
find_orfs.o: find_orfs.cpp find_orfs.h $(HEADERS)

# 3.1
align_seq: align_seq.o $(OBJS)
	$(CPP) -o $@ align_seq.o $(LIBS)
align_seq.o: align_seq.cpp align_seq.h $(HEADERS)

# 3.2
search_seqdb: search_seqdb.o $(OBJS)
	$(CPP) -o $@ search_seqdb.o $(LIBS)
search_seqdb.o: search_seqdb.cpp search_seqdb.h $(HEADERS)

# 4.1
find_motif: find_motif.o $(OBJS)
	$(CPP) -o $@ find_motif.o $(LIBS)
find_motif.o: find_motif.cpp find_motif.h $(HEADERS)

# 4.2
predict_tm: predict_tm.o $(OBJS)
	$(CPP) -o $@ predict_tm.o $(LIBS)
predict_tm.o: predict_tm.cpp predict_tm.h $(HEADERS)

# 5.1
sum_pdb: sum_pdb.o $(OBJS)
	$(CPP) -o $@ sum_pdb.o $(LIBS)
sum_pdb.o: sum_pdb.cpp sum_pdb.h $(HEADERS)

# 5.2
show_pdb_molinfo: show_pdb_molinfo.o $(OBJS)
	$(CPP) -o $@ show_pdb_molinfo.o $(LIBS)
show_pdb_molinfo.o: show_pdb_molinfo.cpp show_pdb_molinfo.h $(HEADERS)

# 6.1
sum_sdf: sum_sdf.o $(OBJS)
	$(CPP) -o $@ sum_sdf.o $(LIBS)
sum_sdf.o: sum_sdf.cpp sum_sdf.h $(HEADERS)

# 6.2
show_atom_types: show_atom_types.o $(OBJS)
	$(CPP) -o $@ show_atom_types.o $(LIBS)
show_atom_types.o: show_atom_types.cpp show_atom_types.h $(HEADERS)

# 7.1
define_secondary_structures: define_secondary_structures.o $(OBJS)
	$(CPP) -o $@ define_secondary_structures.o $(LIBS)
define_secondary_structures.o: define_secondary_structures.cpp\
	define_secondary_structures.h $(HEADERS)

# 7.2
calc_asas: calc_asas.o $(OBJS)
	$(CPP) -o $@ calc_asas.o $(LIBS)
calc_asas.o: calc_asas.cpp calc_asas.h $(HEADERS)

# 8.1
test_superposer: test_superposer.o $(OBJS)
	$(CPP) -o $@ test_superposer.o $(LIBS)
test_superposer.o: test_superposer.cpp test_superposer.h $(HEADERS)

# 8.2
compare_structures: compare_structures.o $(OBJS)
	$(CPP) -o $@ compare_structures.o $(LIBS)
compare_structures.o: compare_structures.cpp compare_structures.h $(HEADERS)

#
# supplementary applications
#
pdb_to_cba: pdb_to_cba.o $(OBJS)
	$(CPP) -o $@ pdb_to_cba.o $(LIBS)
pdb_to_cba.o: pdb_to_cba.cpp pdb_to_cba.h $(HEADERS)

#
count_seqs: count_seqs.o $(OBJS)
	$(CPP) -o $@ count_seqs.o $(LIBS)
count_seqs.o: count_seqs.cpp count_seqs.h $(HEADERS)

#
local_time: local_time.o $(OBJS)
	$(CPP) -o $@ local_time.o $(LIBS)
local_time.o: local_time.cpp local_time.h $(HEADERS)

#
find_identical_mols: find_identical_mols.o $(OBJS)
	$(CPP) -o $@ find_identical_mols.o $(LIBS)
find_identical_mols.o: find_identical_mols.cpp find_identical_mols.h $(HEADERS)

#
morgan_numbering: morgan_numbering.o $(OBJS)
	$(CPP) -o $@ morgan_numbering.o $(LIBS)
morgan_numbering.o: morgan_numbering.cpp morgan_numbering.h $(HEADERS)

#
get_pdb_atom_aas: get_pdb_atom_aas.o $(OBJS)
	$(CPP) -o $@ get_pdb_atom_aas.o $(LIBS)
get_pdb_atom_aas.o: get_pdb_atom_aas.cpp get_pdb_atom_aas.h $(HEADERS)

#
mol_fin_make: mol_fin_make.o $(OBJS)
	$(CPP) -o $@ mol_fin_make.o $(LIBS)
mol_fin_make.o: mol_fin_make.cpp mol_fin_make.h $(HEADERS)

#
mol_fin_cmp: mol_fin_cmp.o $(OBJS)
	$(CPP) -o $@ mol_fin_cmp.o $(LIBS)
mol_fin_cmp.o: mol_fin_cmp.cpp mol_fin_cmp.h $(HEADERS)

#AltPS

altps: altps.o $(OBJS)
	$(CPP) -o $@ altpsm.o $(LIBS)
altps.o: altpsm.cpp altpsm.h $(HEADERS)
#altps.o: altpsm.cpp altpsm.h $(HEADERS) $(LIBCBA) $(LIBMAIN)


######################################################################
# utilities
######################################################################

$(LIBMAIN): $(APP_MAIN)
	$(AR) $(ARFLAGS) $@ $?

$(LIBCBA): $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	-cp -p $(HEADERS) $(INCDIR)

# 1.2
main.o: app.h
app.o: app.cpp app.h

# 2.1
bio_sequence.o: bio_sequence.cpp bio_sequence.h region.h
fasta_reader.o: fasta_reader.cpp fasta_reader.h
bio_seq_cmp.o: bio_seq_cmp.cpp bio_seq_cmp.h bio_sequence.h
bio_seq_cmp_identical.o: bio_seq_cmp_identical.cpp bio_seq_cmp_identical.h\
	bio_seq_cmp.h bio_seq_equiv_region.h bio_sequence.h

# 2.2
dna_sequence.o: dna_sequence.cpp dna_sequence.h bio_sequence.h
genetic_code.o: genetic_code.cpp genetic_code.h

# 3.1
bio_seq_score_matrix.o: bio_seq_score_matrix.cpp bio_seq_score_matrix.h
bio_seq_alignment.o: bio_seq_alignment.cpp bio_seq_alignment.h bio_seq_equiv_region.h
bio_seq_cmp_dp.o: bio_seq_cmp_dp.cpp bio_seq_cmp_dp.h bio_seq_cmp.h\
	bio_sequence.h bio_seq_score_matrix.h bio_seq_alignment.h

# 3.2
bio_seq_db_search.o: bio_seq_db_search.cpp bio_seq_db_search.h bio_sequence.h fasta_reader.h

# 4.2
hydropathy.o: hydropathy.cpp hydropathy.h
trans_membrane_prediction.o: trans_membrane_prediction.cpp trans_membrane_prediction.h\
	bio_sequence.h hydropathy.h

# 5.1
pdb_record.o: pdb_record.cpp pdb_record.h
pdb_entry.o: pdb_entry.cpp pdb_entry.h pdb_record.h molecule.h protein.h
pdb_reader.o: pdb_reader.cpp pdb_reader.h pdb_entry.h pdb_record.h

# 5.2
atom.o: atom.cpp atom.h position.h
position.o: position.cpp position.h matrix.h
matrix.o: matrix.cpp matrix.h
molecule.o: molecule.cpp molecule.h atom.h matrix.h asa_calculator.h
protein.o: protein.cpp protein.h molecule.h atom.h bio_sequence.h

# 6.1
sdf_entry.o: sdf_entry.cpp sdf_entry.h molecule.h
sdf_reader.o: sdf_reader.cpp sdf_reader.h sdf_entry.h

# 7.1
secondary_structure.o: secondary_structure.cpp secondary_structure.h protein.h position.h bit_string.h

# 7.2
asa_calculator.o: asa_calculator.cpp asa_calculator.h atom.h position.h

# 8.1
superposer.o: superposer.cpp superposer.h position.h matrix.h my_math.h

# 8.2
structure_comparison.o: structure_comparison.cpp structure_comparison.h\
	molecule.h protein.h superposer.h

# A1
bit_string.o: bit_string.cpp bit_string.h

# A2
my_math.o: my_math.cpp my_math.h matrix.h

#
cba_entry.o: cba_entry.cpp cba_entry.h atom.h molecule.h protein.h
dir_scan.o: dir_scan.cpp dir_scan.h
pdb_file_processor.o: pdb_file_processor.cpp pdb_file_processor.h dir_scan.h
ncbi_nrpdb.o: ncbi_nrpdb.cpp ncbi_nrpdb.h exception.h
bond_maker.o: bond_maker.cpp bond_maker.h molecule.h protein.h atom.h position.h
mol2_entry.o: mol2_entry.cpp mol2_entry.h molecule.h atom.h
mol2_reader.o: mol2_reader.cpp mol2_reader.h mol2_entry.h

#
date_time.o: date_time.cpp date_time.h

#AltPS
locarea.o: locarea.cpp locarea.h
altpsm.o: altpsm.cpp altpsm.h locarea.h

install:
	-cp -p $(APPS) $(BINDIR)

clean:
	rm -f *.o *.exe $(APPS)
