package simulator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;

public class FileRead {
	
	public static ArrayList<ArrayList<Exon>> loadFiles(String exon_gtf_file, String ref_genome_file){
		ArrayList<ArrayList<Exon>> out = new ArrayList<ArrayList<Exon>>();
		for (int i = 0; i < 25; i++) {
			ArrayList<Exon> temp = new ArrayList<Exon>();
			out.add(temp);
		}
		String temp_line = null;
		Exon exon = null;

		
		BufferedReader exon_gtf = null;
		File read_file = new File(exon_gtf_file);
		try {
			exon_gtf = new BufferedReader(new FileReader(read_file));
			while((temp_line = exon_gtf.readLine()) != null) {
				String[] cols = temp_line.split("\t");
				if (cols[2].equals("exon")) {
					int exon_start = Integer.parseInt(cols[3]);
					int exon_end = Integer.parseInt(cols[4]);
					exon = new Exon();
					exon.setChr_symbol(cols[0]);
					exon.setStart(exon_start);
					exon.setEnd(exon_end);
					int chr_array = exon.getChr_num() - 1;
					if (chr_array >= 0 && chr_array < 25) {
						ArrayList<Exon> exon_list = out.get(chr_array);
						int index = exon.searchMinNoLess(exon_start, exon_list);
						if (index > 0 && exon_start <= exon_list.get(index - 1).getEnd()) {
							Exon front_exon = exon_list.get(index - 1);
							if (exon_end > front_exon.getEnd()) {
								front_exon.setEnd(exon_end);
								while (exon_list.size() > index && exon_end >= exon_list.get(index).getStart()) {
									if (exon_list.get(index).getEnd() > front_exon.getEnd()) {
										front_exon.setEnd(exon_list.get(index).getEnd());
									}
									exon_list.remove(index);
								}
							}
						}
						else if (exon_list.size() > index && exon_end >= exon_list.get(index).getStart()) {
							Exon behind_exon = exon_list.get(index);
							behind_exon.setStart(exon_start);
							if (exon_end > behind_exon.getEnd()) {
								behind_exon.setEnd(exon_end);
								while (exon_list.size() > index + 1 && exon_end >= exon_list.get(index + 1).getStart()) {
									if (exon_list.get(index + 1).getEnd() > behind_exon.getEnd()) {
										behind_exon.setEnd(exon_list.get(index + 1).getEnd());
									}
									exon_list.remove(index + 1);
								}
							}
						}
						else {
							exon_list.add(index, exon);
						}
					}
				}
			}
			exon_gtf.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		BufferedReader genome_fasta = null;
		read_file = new File(ref_genome_file);
		try {
			genome_fasta = new BufferedReader(new FileReader(read_file));
			int chr_array = -1;
			while((temp_line = genome_fasta.readLine()) != null) {
				if (temp_line.charAt(0) == '>') {
					chr_array = exon.chrSymbolToNum(temp_line.substring(1)) - 1;
				}
				else if (chr_array >= 0) {
					ArrayList<Exon> exon_list = out.get(chr_array);
					for (int i = 0; i < exon_list.size(); i++) {
						exon = exon_list.get(i);
						exon.setBase_seq(temp_line.substring(exon.getStart() - 1, exon.getEnd()));
					}
				}
			}
			genome_fasta.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		for (int i = 0; i < out.size();) {
			if (out.get(i).size() > 0) {
				i++;
			}
			else {
				out.remove(i);
			}
		}
		return out;
	}
	
	public static void loadTranscriptsFile(ReadInputs args, ArrayList<Chromosome> out){
		String temp_line = null;
		Chromosome the_chr = null;
		Gene the_gene = null;
		Transcript the_script = null;
		Exon the_exon = null;
		
		BufferedReader exon_gtf = null;
		File read_file = new File(args.getRef_exon_gtf());
		try {
			exon_gtf = new BufferedReader(new FileReader(read_file));
			while((temp_line = exon_gtf.readLine()) != null) {
				String[] cols = temp_line.split("\t");
				if (cols.length < 9) {
					continue;
				}
				if (cols[2].equals("exon")) {
					int exon_start = Integer.parseInt(cols[3]);
					int exon_end = Integer.parseInt(cols[4]);
					the_exon = new Exon();
					the_exon.setChr_symbol(cols[0]);
					the_exon.setStart(exon_start);
					the_exon.setEnd(exon_end);
					the_exon.setScript(the_script);
					if (the_script != null) {
						the_script.addExon(the_exon);
						the_script.setBase_sum(the_script.getBase_sum() + exon_end - exon_start + 1);
					}
				}
				else if (cols[2].equals("transcript")) {
					if (the_script != null) {
						the_script.sortExons(true);
					}
					int script_start = Integer.parseInt(cols[3]);
					int script_end = Integer.parseInt(cols[4]);
					String[] infos = cols[8].split("\"");
					if (infos[2].contains("transcript_id")) {
						if (the_gene != null && the_gene.getId().equals(infos[1])) {
							the_script = new Transcript(the_gene, infos[3], new ArrayList<Exon>(), script_start, script_end, 0, false);
							the_gene.addScript(the_script);
							if (the_chr.getAll_junctions().size() > 0) {
								int index = Method.searchMinNoLess(script_start, the_chr.getAll_junctions());
								if (index == -1) {
									continue;
								}
								int key = 0;
								while(index < the_chr.getAll_junctions().size() && (key = the_chr.getAll_junctions().get(index)) <= script_end) {
									if (index == 0 || key != the_chr.getAll_junctions().get(index - 1)) {
										ArrayList<BackJunction> bj_list = the_chr.getJunction_map().get(key);
										for (int i=0; i < bj_list.size(); i++) {
											if (bj_list.get(i).getStart() >= script_start) {
												bj_list.get(i).addStart_script(the_script);
												if (bj_list.get(i).getEnd() <= script_end) {
													bj_list.get(i).addEnd_script(the_script);
													bj_list.get(i).addBoth_script(the_script);
												}
											}
											else if (bj_list.get(i).getEnd() <= script_end) {
												bj_list.get(i).addEnd_script(the_script);
											}
										}
									}
									index++;
								}
							}
						}
					}
				}
				else if (cols[2].equals("gene")) {
					int gene_start = Integer.parseInt(cols[3]);
					int gene_end = Integer.parseInt(cols[4]);
					String[] infos = cols[8].split("\"");
					if (infos[0].contains("gene_id")) {
						if (the_chr == null || !the_chr.getId().equals(cols[0])) {
							int chr_num = Chromosome.chrSymbolToNum(cols[0]) - 1;
							if (out.get(chr_num) != null) {
								the_chr = out.get(chr_num);
							}
							else{
								the_chr = new Chromosome(cols[0], new ArrayList<Gene>(), new ArrayList<Integer>(), new HashMap<>());
								out.set(chr_num, the_chr);
							}
						}
						the_gene = new Gene(the_chr, infos[1], new ArrayList<Transcript>(), gene_start, gene_end);
						the_chr.addGene(the_gene);
					}
				}
			}
			the_script.sortExons(true);
			exon_gtf.close();
			
			BufferedReader genome_fasta = null;
			read_file = new File(args.getRef_genome_file());
			try {
				genome_fasta = new BufferedReader(new FileReader(read_file));
				int chr_array = -1;
				while((temp_line = genome_fasta.readLine()) != null) {
					if (temp_line.charAt(0) == '>') {
						chr_array = Chromosome.chrSymbolToNum(temp_line.substring(1)) - 1;
					}
					else if (chr_array >= 0) {
						Chromosome chr = out.get(chr_array);
						for (int i = 0; i < chr.getGenes().size(); i++) {
							for (int j = 0; j < chr.getGene(i).getScripts().size(); j++) {
								for (int k = 0; k < chr.getGene(i).getScript(j).getExons().size(); k++) {
									Exon exon = chr.getGene(i).getScript(j).getExon(k);
									exon.setBase_seq(temp_line.substring(exon.getStart() - 1, exon.getEnd()));
								}
							}
						}
						for (Entry<Integer, ArrayList<BackJunction>> bj_list : chr.getJunction_map().entrySet()) {
							for (int i=0; i < bj_list.getValue().size(); i++) {
								BackJunction bj = bj_list.getValue().get(i);
								int length = args.getRead_length();
//								length += length/2;
								bj.setStart_seq(temp_line.substring(bj.getStart() - 1, bj.getStart() + length));
								bj.setEnd_seq(temp_line.substring(bj.getEnd() - length, bj.getEnd() + 1));
							}
						}
					}
				}
				genome_fasta.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			genome_fasta.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static ArrayList<Chromosome> loadJunctionFile(String bed_file){
		ArrayList<Chromosome> out = new ArrayList<>(25);
		for (int i=0; i < 25; i++) {
			out.add(null);
		}
		String temp_line = null;
		BufferedReader reader = null;
		File read_file = new File(bed_file);
		try {
			reader = new BufferedReader(new FileReader(read_file));
			while((temp_line = reader.readLine()) != null) {
				String[] cols = temp_line.split("\t");
				int chr_num = Chromosome.chrSymbolToNum(cols[0]) - 1;
				if (chr_num < 0) {
					continue;
				}
				int start = Integer.parseInt(cols[1]);
				int end = Integer.parseInt(cols[2]);
				Chromosome chr = null;
				if (out.get(chr_num) == null) {
					chr = new Chromosome(cols[0], new ArrayList<Gene>(), new ArrayList<Integer>(), new HashMap<Integer, ArrayList<BackJunction>>());
					out.set(chr_num, chr);
				}
				else {
					chr = out.get(chr_num);
				}
				BackJunction bj = new BackJunction();
				bj.setStart(start);
				bj.setEnd(end);
				chr.getAll_junctions().add(start);
				chr.getAll_junctions().add(end);
				if (chr.getJunction_map().containsKey(start)) {
					chr.getJunction_map().get(start).add(bj);
				}
				else {
					ArrayList<BackJunction> bj_list = new ArrayList<>();
					bj_list.add(bj);
					chr.getJunction_map().put(start, bj_list);
				}
				if (chr.getJunction_map().containsKey(end)) {
					chr.getJunction_map().get(end).add(bj);
				}
				else {
					ArrayList<BackJunction> bj_list = new ArrayList<>();
					bj_list.add(bj);
					chr.getJunction_map().put(end, bj_list);
				}
			}
			for (int i=0; i < 25; i++) {
				if (out.get(i) != null) {
					Collections.sort(out.get(i).getAll_junctions());
				}
			}

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return out;
	}
}
