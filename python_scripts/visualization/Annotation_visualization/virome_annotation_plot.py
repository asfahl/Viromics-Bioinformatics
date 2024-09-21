from pycirclize import Circos
from pycirclize.parser import Gff
import os, sys

def main():

	prodigal_filename = os.path.abspath(sys.argv[1])
	assert prodigal_filename.endswith(".gff")
	
	phanotate_filename = os.path.abspath(sys.argv[2])
	assert phanotate_filename.endswith(".gff")
	
	out_filename = os.path.abspath(sys.argv[3])
	assert out_filename.endswith(".png")
	
	contig_id = sys.argv[4]
	
	prodigal_gff = Gff(prodigal_filename)
	phanotate_gff = Gff(phanotate_filename)

	# Initialize circos sector by genome size
	circos = Circos(sectors={phanotate_gff.name: phanotate_gff.range_size})
	circos.text(contig_id, size=15)
	sector = circos.sectors[0]
	
	# Outer track
	outer_track = sector.add_track((98, 100))
	outer_track.axis(fc="lightgrey")
	outer_track.xticks_by_interval(5000, label_formatter=lambda v: f"{v / 1000:.0f} Kb")
	outer_track.xticks_by_interval(1000, tick_length=1, show_label=False)
	
	# Plot forward & reverse CDS genomic features
	cds_track = sector.add_track((90, 95))
	cds_track.genomic_features(phanotate_gff.extract_features("CDS", target_strand=1), plotstyle="arrow", fc="salmon")
	cds_track.genomic_features(phanotate_gff.extract_features("CDS", target_strand=-1), plotstyle="arrow", fc="skyblue")
	
	cds_track = sector.add_track((82, 87))
	cds_track.genomic_features(prodigal_gff.extract_features("CDS", target_strand=1), plotstyle="arrow", fc="salmon")
	cds_track.genomic_features(prodigal_gff.extract_features("CDS", target_strand=-1), plotstyle="arrow", fc="skyblue")
	
	circos.savefig(out_filename)
	
if __name__ == "__main__":
    main()

