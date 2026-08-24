##
## Makefile for all executables
##

## Default compilation flags.
## Override with:
##   make CXXFLAGS=XXXXX
## -include cstdint makes fixed-width integer types (int32_t, uint32_t, ...) available
## on strict modern compilers (e.g. gcc 13+), which no longer pull them in transitively.
CXXFLAGS= -O3 -g -D__STDC_LIMIT_MACROS -D_FILE_OFFSET_BITS=64 -std=c++0x -DMACOSX -pthread -include cstdint

## To create a static distribution file, run:
##   make static-dist
ifeq ($(STATIC),1)
LDFLAGS=-static
else
LDFLAGS=
endif

## Portable distribution build (opt-in). Produces a fully-static binary with no
## network dependencies (no libcurl/OpenSSL) and portable SSE4.1 SIMD for spoa, so the
## result runs on any x86-64 Linux. Orthogonal to STATIC above: a plain `make` and
## `make STATIC=1` are byte-for-byte unaffected. Build with:  make PORTABLE=1
NETWORK_LIBS        = -lcurl -lcrypto
HTSLIB_CONFIG_FLAGS = --enable-gcs --enable-s3 --enable-libcurl
SPOA_CMAKE_FLAGS    =
ifeq ($(PORTABLE),1)
LDFLAGS             = -static
NETWORK_LIBS        =
HTSLIB_CONFIG_FLAGS = --disable-libcurl --disable-gcs --disable-s3 --disable-plugins
SPOA_CMAKE_FLAGS    = -Dspoa_optimize_for_native=OFF -Dspoa_optimize_for_portability=ON
endif



## Source code files, add new files to this list
SRC_COMMON  = src/base_quality.cpp src/error.cpp src/region.cpp src/stringops.cpp src/zalgorithm.cpp src/alignment_filters.cpp src/extract_indels.cpp src/mathops.cpp src/pcr_duplicates.cpp src/bam_io.cpp
SRC_HIPSTR  = src/hipstr_main.cpp src/bam_processor.cpp src/stutter_model.cpp src/snp_phasing_quality.cpp src/snp_tree.cpp src/em_stutter_genotyper.cpp src/seq_stutter_genotyper.cpp src/snp_bam_processor.cpp src/genotyper_bam_processor.cpp src/vcf_input.cpp src/read_pooler.cpp src/version.cpp src/haplotype_tracker.cpp src/pedigree.cpp src/vcf_reader.cpp src/genotyper.cpp src/directed_graph.cpp src/debruijn_graph.cpp src/fasta_reader.cpp src/vcf_writer.cpp
SRC_SEQALN  = src/SeqAlignment/HapAligner.cpp src/SeqAlignment/AlignmentOps.cpp src/SeqAlignment/HapBlock.cpp src/SeqAlignment/NeedlemanWunsch.cpp src/SeqAlignment/Haplotype.cpp src/SeqAlignment/HaplotypeGenerator.cpp src/SeqAlignment/HTMLCreator.cpp src/SeqAlignment/AlignmentViz.cpp src/SeqAlignment/AlignmentTraceback.cpp src/SeqAlignment/StutterAlignerClass.cpp
SRC_DENOVO  = src/denovos/denovo_main.cpp src/error.cpp src/stringops.cpp src/version.cpp src/pedigree.cpp src/haplotype_tracker.cpp src/vcf_input.cpp src/denovos/denovo_scanner.cpp src/mathops.cpp src/vcf_reader.cpp src/denovos/denovo_allele_priors.cpp src/denovos/trio_denovo_scanner.cpp

CEPHES_ROOT=lib/cephes

LIBS              = -L./ -lm -Llib/htslib/lib -lhts -ldeflate -lz -L$(CEPHES_ROOT)/ -llzma -lbz2 $(NETWORK_LIBS) -Llib/spoa/build/lib -lspoa
INCLUDE           = -Ilib -Ilib/htslib/include -Ilib/spoa/include
CEPHES_LIB        = lib/cephes/libprob.a
HTSLIB_LIB        = lib/htslib/lib/libhts.a

# For each CPP file, generate an object file
OBJ_COMMON  := $(SRC_COMMON:.cpp=.o)
OBJ_HIPSTR  := $(SRC_HIPSTR:.cpp=.o)
OBJ_SEQALN  := $(SRC_SEQALN:.cpp=.o)
OBJ_DENOVO  := $(SRC_DENOVO:.cpp=.o)

.PHONY: all
all: HTSLIB-docker SPOA-docker LongTR DenovoFinder test/fast_ops_test test/haplotype_test test/read_vcf_alleles_test test/snp_tree_test test/vcf_snp_tree_test

# Rebuild every object when any header changes. The `include $(subst .cpp,.d,$(SRC))` further down
# keys off an undefined SRC, so it expands to nothing and header dependencies have never been
# tracked: editing a header rebuilt nothing, and the resulting binary could link objects compiled
# against two different layouts of the same class (observed as heap corruption at exit).
# Must stay BELOW `all`, or the first object becomes make's default goal.
$(OBJ_COMMON) $(OBJ_HIPSTR) $(OBJ_SEQALN) $(OBJ_DENOVO): $(wildcard src/*.h src/SeqAlignment/*.h src/denovos/*.h)

# Create a tarball with static binaries
.PHONY: static-dist
static-dist:
	rm -f LongTR
	$(MAKE) STATIC=1
	( VER="$$(git describe --abbrev=7 --dirty --always --tags)" ;\
	  DST="LongTR-$${VER}-static-$$(uname -s)-$$(uname -m)" ; \
	  mkdir "$${DST}" && \
            mkdir "$${DST}/scripts" && \
            cp LongTR VizAln VizAlnPdf README.md "$${DST}" && \
            cp scripts/filter_haploid_vcf.py scripts/filter_vcf.py scripts/generate_aln_html.py scripts/html_alns_to_pdf.py "$${DST}/scripts" && \
            tar -czvf "$${DST}.tar.gz" "$${DST}" && \
            rm -r "$${DST}/" \
        )

# Create a tarball with a portable, fully-static binary (opt-in; no network deps,
# portable SIMD). Intended for distributable releases. Usage:  make portable-dist
.PHONY: portable-dist
portable-dist:
	rm -f LongTR
	$(MAKE) PORTABLE=1 HTSLIB-docker SPOA-docker LongTR DenovoFinder && strip LongTR DenovoFinder
	( VER="$${VER:-$$(git describe --abbrev=7 --dirty --always --tags 2>/dev/null || echo dev)}" ;\
	  DST="LongTR-$${VER}-portable-$$(uname -s)-$$(uname -m)" ; \
	  mkdir "$${DST}" && \
            mkdir "$${DST}/scripts" && \
            cp LongTR VizAln VizAlnPdf README.md "$${DST}" && \
            cp scripts/filter_haploid_vcf.py scripts/filter_vcf.py scripts/generate_aln_html.py scripts/html_alns_to_pdf.py "$${DST}/scripts" && \
            tar -czvf "$${DST}.tar.gz" "$${DST}" && \
            rm -r "$${DST}/" \
        )

version:
	git describe --abbrev=7 --dirty --always --tags | awk '{print "#include \"version.h\""; print "const std::string VERSION = \""$$0"\";"}' > src/version.cpp

# Clean the generated files of the main project only
.PHONY: clean
clean:
	rm -f *~ src/*.o src/*.d src/*~ src/SeqAlignment/*~ src/SeqAlignment/*.o src/denovos/*~ src/denovos/*.o LongTR DenovoFinder test/allele_expansion_test test/fast_ops_test test/haplotype_test test/read_vcf_alleles_test test/snp_tree_test test/vcf_snp_tree_test

# Clean all compiled files
.PHONY: clean-all
clean-all: clean
	cd lib/htslib && $(MAKE) clean
	rm lib/cephes/*.o $(CEPHES_LIB)

# The GNU Make trick to include the ".d" (dependencies) files.
# If the files don't exist, they will be re-generated, then included.
# If this causes problems with non-gnu make (e.g. on MacOS/FreeBSD), remove it.
include $(subst .cpp,.d,$(SRC))

# The resulting binary executable

.PHONY: HTSLIB
HTSLIB:
	@if [ ! -d "lib/htslib" ]; then \
		cd lib && git clone --recurse-submodules https://github.com/samtools/htslib.git && cd ..;\
	else\
		echo "htslib directory already exists in lib/ folder";\
	fi
.PHONY: HTSLIB-update
HTSLIB-update: HTSLIB 
	@cd lib/htslib && git pull

.PHONY: HTSLIB-docker
HTSLIB-docker: HTSLIB-update
	@cd lib/htslib && autoreconf -i && ./configure --prefix="$(CURDIR)"/lib/htslib $(HTSLIB_CONFIG_FLAGS) && make -j && make install

.PHONY: SPOA
SPOA:
	@if [ ! -d "lib/spoa" ]; then \
		cd lib && git clone https://github.com/rvaser/spoa.git && cd ..;\
	else\
		echo "spoa directory already exists in lib/ folder";\
	fi

.PHONY: SPOA-update
SPOA-update: SPOA
	@cd lib/spoa && git pull

.PHONY: SPOA-docker
SPOA-docker: SPOA-update
	@if [ ! -d "lib/spoa/build" ]; then \
		cd lib/spoa && mkdir build && cd build && cmake $(SPOA_CMAKE_FLAGS) -DCMAKE_BUILD_TYPE=Release .. && cd .. && make -C build;\
	else\
		cd lib/spoa/build && cmake $(SPOA_CMAKE_FLAGS) -DCMAKE_BUILD_TYPE=Release .. && cd .. && make -C build;\
	fi

LongTR: $(OBJ_COMMON) $(OBJ_HIPSTR) $(HTSLIB_LIB) $(OBJ_SEQALN)
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

DenovoFinder: $(OBJ_DENOVO) $(HTSLIB_LIB)
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

PhasingChecker: src/check_phasing.cpp src/region.cpp src/error.cpp src/haplotype_tracker.cpp src/version.cpp src/pedigree.cpp src/vcf_reader.cpp src/stringops.cpp $(HTSLIB_LIB)
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

test/haplotype_test: test/haplotype_test.cpp src/SeqAlignment/Haplotype.cpp src/SeqAlignment/HapBlock.cpp src/SeqAlignment/NeedlemanWunsch.cpp src/error.cpp src/stringops.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

test/em_stutter_test: test/em_stutter_test.cpp src/em_stutter_genotyper.cpp src/genotyper_bam_processor.cpp src/error.cpp src/mathops.cpp src/stringops.cpp src/stutter_model.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

test/fast_ops_test: test/fast_ops_test.cpp src/mathops.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^

test/read_vcf_alleles_test: test/read_vcf_alleles_test.cpp src/error.cpp src/region.cpp src/vcf_input.cpp src/vcf_reader.cpp $(HTSLIB_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

test/snp_tree_test: src/snp_tree.cpp src/error.cpp test/snp_tree_test.cpp src/haplotype_tracker.cpp src/vcf_reader.cpp $(HTSLIB_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

test/vcf_snp_tree_test: test/vcf_snp_tree_test.cpp src/error.cpp src/snp_tree.cpp src/haplotype_tracker.cpp src/vcf_reader.cpp $(HTSLIB_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

# Build each object file independently
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

# Auto-Generate header dependencies for each CPP file.
%.d: %.cpp
	$(CXX) -c -MP -MD $(CXXFLAGS) $(INCLUDE) $< > $@

# Rebuild CEPHES library if needed
$(CEPHES_LIB):
	cd lib/cephes && $(MAKE) -fPIE -pie

