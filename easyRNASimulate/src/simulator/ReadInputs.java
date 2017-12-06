package simulator;

public class ReadInputs {
	
	private String ref_genome_file=null;
	private String ref_exon_gtf=null;
	private String out_file=null;
	private String out_dir=null;
	private int sim_mode=0;
	private int read_count=0;
	private int read_length=0;
	private int min_splice=10;
	private int max_splice=100;
	private double circ_scale=0.1;
	private double splice_input_scale=0.1;
	private boolean exon_bound=true;
	private boolean read_pair=false;
	private int FULL_NOTE=2047;
	
	public ReadInputs() {
	}
	
	public ReadInputs(String[] args) {
		for (int i = 0; i < args.length; i++) {
			if (args[i].charAt(0) == '-') {
				String temp_args = args[i].toLowerCase();
				if (temp_args.equals("-g")) {
					i++;
					this.ref_genome_file=args[i];
				}
				else if(temp_args.equals("-e")) {
					i++;
					this.ref_exon_gtf=args[i];
				}
				else if(temp_args.equals("-o")) {
					i++;
					this.out_file=args[i];
				}
				else if(temp_args.equals("-m")) {
					i++;
					this.sim_mode=Integer.parseInt(args[i]);
				}
				else if(temp_args.equals("-c")) {
					i++;
					this.read_count=Integer.parseInt(args[i]);
				}
				else if(temp_args.equals("-l")) {
					i++;
					this.read_length=Integer.parseInt(args[i]);
				}
				else if(temp_args.equals("-mins")) {
					i++;
					this.min_splice=Integer.parseInt(args[i]);
				}
				else if(temp_args.equals("-maxs")) {
					i++;
					this.max_splice=Integer.parseInt(args[i]);
				}
				else if(temp_args.equals("-circ")) {
					i++;
					this.circ_scale=Double.parseDouble(args[i]);
				}
				else if(temp_args.equals("-s")) {
					i++;
					this.splice_input_scale=Double.parseDouble(args[i]);
				}
				else if(temp_args.equals("-b")) {
					this.exon_bound=false;
				}
				else if(temp_args.equals("-p")) {
					this.read_pair=true;
				}
				else if(temp_args.equals("-h")) {
					this.printNote(this.FULL_NOTE);
				}
				else {
					System.out.println("Unknown parameter: " + args[i]);
				}
			}
			else {
				System.out.println("Unknown parameter: " + args[i]);
			}
		}
	}
	
	int completeArgs() {
		int out=0;
		if (this.ref_genome_file == null) {
			out |= 1;
		}
		if (this.ref_exon_gtf == null) {
			out |= 2;
		}
		if (this.out_file == null && this.out_dir == null) {
			out |= 4;
		}
		if (this.sim_mode <= 0 || this.sim_mode > SimulateModes.getTotal_modes()) {
			out |= 8;
		}
		if (this.read_count <= 0) {
			out |= 16;
		}
		if (this.read_length <= 0) {
			out |= 32;
		}
		if (this.circ_scale < 0.0 || this.circ_scale > 1.0) {
			out |= 64;
		}
		if (this.splice_input_scale < 0.0 || this.splice_input_scale > 1.0) {
			out |= 128;
		}
		return out;
	}
	
	void printNote(int note) {
		String[] note_lines= {
				"\t-g\trequires ref genome fasta file",
				"\t-e\trequires ref exon gtf file",
				"\t-o/-d\trequires output filename/directory",
				"\t-m\trequires integer for choosing mode",
				"\t-c\trequires integer for amount of reads",
				"\t-l\trequires integer for read length",
				"\t-circ\tneeds double for scale of circ reads in spliced reads from 0.0 to 1.0",
				"\t-s\tneeds double for scale of spliced reads from 0.0 to 1.0",
				"\t-b\tdisable splice from the exon boundary",
				"\t-p\tenable pair end reads",
				"\t-h\tdisplay the help note"
		};
		if (note == this.FULL_NOTE) {
			System.out.println("This is help document");
		}
		for (int i = 0; i < note_lines.length; i++) {
			if ((note >>> i & 1) == 1) {
				System.out.println(note_lines[i]);
			}
		}
	}
	
	String getRef_genome_file() {
		return ref_genome_file;
	}
	
	String getRef_exon_gtf() {
		return ref_exon_gtf;
	}
	
	String getOut_file() {
		return out_file;
	}
	
	String getOut_dir() {
		return out_dir;
	}
	
	int getSim_mode() {
		return sim_mode;
	}
	
	int getRead_count() {
		return read_count;
	}
	
	int getRead_length() {
		return read_length;
	}
	
	int getMin_splice() {
		return min_splice;
	}

	int getMax_splice() {
		return max_splice;
	}

	double getCirc_scale() {
		return circ_scale;
	}
	
	double getSplice_input_scale() {
		return splice_input_scale;
	}
	
	boolean isExon_bound() {
		return exon_bound;
	}

	boolean isRead_pair() {
		return read_pair;
	}

	int getFULL_NOTE() {
		return FULL_NOTE;
	}
	
}
