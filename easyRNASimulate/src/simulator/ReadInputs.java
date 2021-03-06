package simulator;

public class ReadInputs {
	
	private String ref_genome_file=null;
	private String ref_exon_gtf=null;
	private String ref_bed_file=null;
	private String out_file=null;
	private String out_dir=null;
	private int sim_mode=5;
	private int read_count=0;
	private int read_length=0;
//	private int min_splice=-1;
//	private int max_splice=-1;
	private int peak_num = 2000;
	private int circ_num = 20000;
//	private double circ_scale=0.1;
//	private double splice_input_scale=0.1;
	private double depth = 30.0;
	private double enrich = 100.0;
	private double ip_scale = 0.5;
	private boolean exon_bound=true;
	private boolean read_pair=false;
	private int FULL_NOTE=65535;

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
				else if(temp_args.equals("-bed")) {
					i++;
					this.ref_bed_file=args[i];
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
//				else if(temp_args.equals("-mins")) {
//					i++;
//					this.min_splice=Integer.parseInt(args[i]);
//				}
//				else if(temp_args.equals("-maxs")) {
//					i++;
//					this.max_splice=Integer.parseInt(args[i]);
//				}
				else if(temp_args.equals("-peak")) {
					i++;
					this.peak_num=Integer.parseInt(args[i]);
				}
				else if(temp_args.equals("-circ")) {
					i++;
					this.circ_num=Integer.parseInt(args[i]);
				}
//				else if(temp_args.equals("-s")) {
//					i++;
//					this.splice_input_scale=Double.parseDouble(args[i]);
//				}
				else if(temp_args.equals("-dep")) {
					i++;
					this.depth=Double.parseDouble(args[i]);
				}
				else if(temp_args.equals("-enrich")) {
					i++;
					this.enrich=Double.parseDouble(args[i]);
				}
				else if(temp_args.equals("-is")) {
					i++;
					this.ip_scale=Double.parseDouble(args[i]);
				}
				else if(temp_args.equals("-no_b")) {
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
//		if (this.min_splice < 0) {
//			this.min_splice = (int) (this.read_length * 0.1);
//		}
//		if (this.max_splice < 0) {
//			this.max_splice = (int) (this.read_length * 0.8);
//		}
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
//		if (this.circ_scale < 0.0 || this.circ_scale > 1.0) {
//			out |= 64;
//		}
//		if (this.splice_input_scale < 0.0 || this.splice_input_scale > 1.0) {
//			out |= 128;
//		}
//		if (this.min_splice > this.read_length) {
//			out |= 256;
//		}
//		if (this.max_splice > this.read_length) {
//			out |= 512;
//		}
		return out;
	}
	
	void printNote(int note) {
		String[] warn_lines= {
				"\t-g\trequires ref genome fasta file",
				"\t-e\trequires ref exon gtf file",
				"\t-o/-d\trequires output filename/directory",
				"\t-m\trequires integer for choosing mode",
				"\t-c\trequires integer for amount of reads",
				"\t-l\trequires integer for read length",
//				"\t-s\tneeds double for scale of spliced reads from 0.0 to 1.0",
//				"\t-mins\tneeds int for minimum splice segment length (default is 0.1*read length, should less than read length)",
//				"\t-maxs\tneeds int for maximum splice segment length (default is 0.8*read length, should less than read length)",
				
		};

		String help_line = "\t-bed\tneeds circRNA bed file\n"
				+ "\t-peak\tneeds int for total peak count\n"
				+ "\t-circ\tneeds double for scale of circ reads in spliced reads from 0.0 to 1.0\n"
				+ "\t-b\tdisable splice from the exon boundary\n"
				+ "\t-p\tenable pair end reads\n"
				+ "\t-dep\tdouble for read depth (deafault is 30.0)\n"
				+ "\t-enrich\tdouble for peak enrichment (default is 100.0)\n"
				+ "\t-is\tdouble for IP read depth compared with INPUT (default is 0.5)\n"
				+ "\t-h\tdisplay the help note";
		if (note == this.FULL_NOTE) {
			System.out.println("This is help document");
			for (int i = 0; i < warn_lines.length; i++) {
				System.out.println(warn_lines[i]);
			}
			System.out.println(help_line);
		}
		else {
			for (int i = 0; i < warn_lines.length; i++) {
				if ((note >>> i & 1) == 1) {
					System.out.println(warn_lines[i]);
				}
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
	
	String getRef_bed_file() {
		return ref_bed_file;
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
	
//	int getMin_splice() {
//		return min_splice;
//	}
//
//	int getMax_splice() {
//		return max_splice;
//	}

	int getPeak_num() {
		return peak_num;
	}

	public int getCirc_num() {
		return circ_num;
	}

//	double getCirc_scale() {
//		return circ_scale;
//	}
//	
//	double getSplice_input_scale() {
//		return splice_input_scale;
//	}
	
	public double getDepth() {
		return depth;
	}

	public double getEnrich() {
		return enrich;
	}

	public double getIp_scale() {
		return ip_scale;
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
	
	void setRead_count(int read_count) {
		this.read_count = read_count;
	}
	
    protected ReadInputs clone() throws CloneNotSupportedException {  
        return (ReadInputs)super.clone();  
    }  
	
}
