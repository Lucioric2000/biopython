#!/usr/bin/env bash
#
# Usage:
#   ./test.sh [--output_path <junit.xml>] <base|new>
#
#   base  run the existing repository tests before the changes; these
#         must pass both before and after the changes are applied.
#   new   run the new tests; these fail before the changes and pass
#         after them.
set -uo pipefail

cd /app

JUNIT_XML_ARGUMENT=()
MODE=""

while [ "$#" -gt 0 ]; do
  case "$1" in
    --output_path)
      if [ "$#" -lt 2 ]; then
        echo "--output_path requires a value" >&2
        exit 2
      fi

      OUTPUT_PATH="$2"
      if [[ "$OUTPUT_PATH" != /* ]]; then
        OUTPUT_PATH="/app/$OUTPUT_PATH"
      fi

      mkdir -p "$(dirname "$OUTPUT_PATH")"
      JUNIT_XML_ARGUMENT=(--junitxml "$OUTPUT_PATH")
      shift 2
      ;;
    base|new)
      MODE="$1"
      shift
      ;;
    *)
      echo "unknown argument: $1" >&2
      exit 2
      ;;
  esac
done


echo Running $MODE tests with ${JUNIT_XML_ARGUMENT[@]}

#We execute the tests not requiring a real internet connection (like requires test_send_request)
case "$MODE" in
  base)
    cd Tests; python -m pytest test_Ace.py test_Affy.py test_AlignIO.py test_AlignIO_ClustalIO.py test_AlignIO_EmbossIO.py test_AlignIO_FastaIO.py test_AlignIO_MauveIO.py test_AlignIO_PhylipIO.py test_AlignIO_convert.py test_Align_Alignment.py test_Align_a2m.py test_Align_aligncore.py test_Align_bed.py test_Align_bigbed.py test_Align_bigmaf.py test_Align_bigpsl.py test_Align_chain.py test_Align_clustal.py test_Align_codonalign.py test_Align_emboss.py test_Align_exonerate.py test_Align_fasta.py test_Align_format_matrix.py test_Align_hhr.py test_Align_maf.py test_Align_mauve.py test_Align_msf.py test_Align_nexus.py test_Align_phylip.py test_Align_psl.py test_Align_sam.py test_Align_stockholm.py test_Align_tabular.py test_Blast_Record.py test_Blast_parser.py test_CAPS.py test_Chi2.py test_Cluster.py test_CodonTable.py test_ColorSpiral.py test_Compass.py test_Consensus.py test_EMBL_unittest.py test_EmbossPrimer.py test_Enzyme.py test_ExPASy_offline.py test_File.py test_GenBank.py test_GenomeDiagram.py test_GraphicsChromosome.py test_GraphicsDistribution.py test_GraphicsGeneral.py test_KEGG.py test_KGML_graphics.py test_KGML_nographics.py test_KeyWList.py test_MafIO_index.py test_Medline.py test_NMR.py test_Nexus.py test_PDB_CEAligner.py test_PDB_DSSP.py test_PDB_Dice.py test_PDB_Disordered.py test_PDB_Exposure.py test_PDB_FragmentMapper.py test_PDB_KDTree.py test_PDB_MMCIF2Dict.py test_PDB_MMCIFIO.py test_PDB_MMCIFParser.py test_PDB_NACCESS.py test_PDB_PDBIO.py test_PDB_PDBMLParser.py test_PDB_PDBParser.py test_PDB_Polypeptide.py test_PDB_QCPSuperimposer.py test_PDB_ResidueDepth.py test_PDB_SASA.py test_PDB_SMCRA.py test_PDB_Selection.py test_PDB_StructureAlignment.py test_PDB_Superimposer.py test_PDB_binary_cif.py test_PDB_internal_coords.py test_PDB_parse_pdb_header.py test_PDB_vectors.py test_PQR.py test_Pathway.py test_Phd.py test_Phylo.py test_PhyloXML.py test_Phylo_CDAO.py test_Phylo_NeXML.py test_Phylo_igraph.py test_Phylo_matplotlib.py test_Phylo_networkx.py test_PopGen_GenePop_nodepend.py test_ProtParam.py test_RCSBFormats.py test_Restriction.py test_SCOP_Astral.py test_SCOP_Cla.py test_SCOP_Des.py test_SCOP_Dom.py test_SCOP_Hie.py test_SCOP_Raf.py test_SCOP_Residues.py test_SCOP_Scop.py test_SVDSuperimposer.py test_SearchIO_blast_tab.py test_SearchIO_blast_tab_index.py test_SearchIO_blast_xml.py test_SearchIO_blast_xml_index.py test_SearchIO_blat_psl.py test_SearchIO_blat_psl_index.py test_SearchIO_exonerate.py test_SearchIO_exonerate_text_index.py test_SearchIO_exonerate_vulgar_index.py test_SearchIO_fasta_m10.py test_SearchIO_fasta_m10_index.py test_SearchIO_hhsuite2_text.py test_SearchIO_hmmer2_text.py test_SearchIO_hmmer2_text_index.py test_SearchIO_hmmer3_domtab.py test_SearchIO_hmmer3_domtab_index.py test_SearchIO_hmmer3_tab.py test_SearchIO_hmmer3_tab_index.py test_SearchIO_hmmer3_text.py test_SearchIO_hmmer3_text_index.py test_SearchIO_infernal_tab.py test_SearchIO_infernal_tab_index.py test_SearchIO_infernal_text.py test_SearchIO_infernal_text_index.py test_SearchIO_interproscan_xml.py test_SearchIO_model.py test_SearchIO_write.py test_SeqFeature.py test_SeqIO.py test_SeqIO_AbiIO.py test_SeqIO_FastaIO.py test_SeqIO_Gck.py test_SeqIO_Gfa.py test_SeqIO_Insdc.py test_SeqIO_NibIO.py test_SeqIO_PdbIO.py test_SeqIO_QualityIO.py test_SeqIO_SeqXML.py test_SeqIO_SffIO.py test_SeqIO_SnapGene.py test_SeqIO_TwoBitIO.py test_SeqIO_UniprotIO.py test_SeqIO_Xdna.py test_SeqIO_features.py test_SeqIO_index.py test_SeqIO_write.py test_SeqRecord.py test_SeqUtils.py test_Seq_objs.py test_SwissProt.py test_TreeConstruction.py test_Tutorial.py test_UniGene.py test_align.py test_align_substitution_matrices.py test_bgzf.py test_codonalign.py test_geo.py test_mmtf.py test_motifs.py test_pairwise2.py test_pairwise2_no_C.py test_pairwise_aligner.py test_pairwise_alignment_map.py test_phenotype.py test_phenotype_fit.py test_prodoc.py test_prosite.py test_seq.py test_translate.py -k "not test_random and not test_random_sequences and not test_correct_retrieval_1" -v "${JUNIT_XML_ARGUMENT[@]}"
    ;;
  new)
    cd Tests; python -m pytest -v test_PDB_gromacs.py "${JUNIT_XML_ARGUMENT[@]}"
    ;;
  *)
    echo "unknown mode: $MODE (expected base or new)" >&2
    exit 2
    ;;
esac

echo listing the output path $OUTPUT_PATH
ls $OUTPUT_PATH