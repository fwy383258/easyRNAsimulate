package simulator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;


public class FileWrite {

	public void fileWrite(String fileName, ArrayList<String> out){
		File fo = new File (fileName);
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(fo));
			for (String tempString : out){
				writer.write(tempString);
				writer.newLine();
			}
			writer.flush();
			writer.close();
		}
		catch(IOException e) {
			e.printStackTrace();
		}
		finally {
			try {
				writer.close();
			}
			catch(IOException e1){
			}
		}
	}
	
	public void bedWrite(ArrayList<Bed12> out, String fileName){
		File fo = new File (fileName);
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(fo));
			for (Bed12 record : out){
				writer.write(record.toString());
				writer.newLine();
			}
			writer.flush();
			writer.close();
		}
		catch(IOException e) {
			e.printStackTrace();
		}
		finally {
			try {
				writer.close();
			}
			catch(IOException e1){
			}
		}
	}
	
	public void fileAppend(String fileName, ArrayList<String> out){
		File fo = new File (fileName);
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(fo, true));
			for (String tempString : out){
				writer.write(tempString);
				writer.newLine();
			}
			writer.flush();
			writer.close();
		}
		catch(IOException e) {
			e.printStackTrace();
		}
		finally {
			try {
				writer.close();
			}
			catch(IOException e1){
			}
		}
	}
	
	public void exonWrite(String file_name, ArrayList<ArrayList<Exon>> out) {
		File file_out = new File (file_name);
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(file_out));
			for (int i = 0; i < out.size(); i++) {
				int last_end = 0;
				for (int j = 0; j < out.get(i).size(); j++) {
					Exon exon = out.get(i).get(j);
					if (last_end >= exon.getStart()) {
						writer.write("DUP!\t");
					}
					writer.write(exon.getChr_symbol());
					writer.write("\t");
					writer.write(String.valueOf(exon.getStart()));
					writer.write("\t");
					last_end = exon.getEnd();
					writer.write(String.valueOf(last_end));
					writer.write("\t");
					writer.write(exon.getBase_seq());
					writer.newLine();
				}
			}
			writer.flush();
			writer.close();
		}
		catch(IOException e) {
			e.printStackTrace();
		}
		finally {
			try {
				writer.close();
			}
			catch(IOException e) {
			}
		}
	}
	
	public void readsWrite(ReadInputs args, String out_file, ArrayList<Read> chr_read, int chr) {
		BufferedReader reader = null;
		File read_file = new File(args.getRef_genome_file());
		BufferedWriter writer = null;
		String line = null;
		StringBuffer qb = new StringBuffer();
		for (int i = 0; i < args.getRead_length(); i ++) {
			qb.append('I');
		}
		String quality = qb.toString();
		try {
			reader = new BufferedReader(new FileReader(read_file));
			writer = new BufferedWriter(new FileWriter(out_file));
			int chr_array = -1;
			boolean first_line = false;
			while((line = reader.readLine()) != null) {
				if (line.charAt(0) == '>') {
					chr_array = Chromosome.chrSymbolToNum(line.substring(1)) - 1;
				}
				else if (chr_array == chr) {
					if (chr_read == null) {
						continue;
					}
					int read_num=0;
					int chr_length = line.length();
					for (int i = 0; i < chr_read.size(); i++) {
						read_num++;
						Read read = chr_read.get(i);
						StringBuffer id = new StringBuffer();
						StringBuffer base_seq = new StringBuffer();
						id.append('@');
						if (read.isCirc_flag()) {
							id.append("circ:");
						}
						if (read.isPeak_flag()) {
							id.append("peak:");
						}
						id.append(read.getId());
						id.append(':');
						id.append(read_num);
						id.append(':');
						id.append(Chromosome.chrNumToSymbol(chr_array + 1));
						id.append(':');
						for (int j = 0; j < read.getPositions().size(); j += 2) {
							id.append(read.getPositions().get(j));
							id.append('-');
							id.append(read.getPositions().get(j + 1));
							id.append('_');
							base_seq.append(line.substring(read.getPositions().get(j) - 1, read.getPositions().get(j + 1)));
						}
						id.setLength(id.length() - 1);
						int base_end = read.getPositions().get(read.getPositions().size() - 1) + 10;
						if (base_end > chr_length) {
							base_end = chr_length;
						}
						base_seq.append(line.substring(read.getPositions().get(read.getPositions().size() - 1), base_end));
						Method.randomErrorSeq(base_seq);
						if (base_seq.length() >= args.getRead_length()) {
							base_seq.setLength(args.getRead_length());
						}
						else {
							System.out.println("Seq too short");
							int to_add = args.getRead_length() - base_seq.length();
							for (int j = 0; j < to_add; j++) {
								base_seq.append('N');
							}
						}
						if (first_line) {
							writer.newLine();
						}
						writer.write(id.toString());
						writer.newLine();
						writer.write(base_seq.toString());
						writer.newLine();
						writer.write('+');
						writer.newLine();
						writer.write(quality);
						first_line = true;
					}
					writer.flush();
				}
			}
			writer.flush();
			reader.close();
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e1) {
				}
			}
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e1) {
				}
			}
		}
	}
	
	public void readsWrite(ReadInputs args, String out_file, ArrayList<ArrayList<Read>> reads) {
		BufferedReader reader = null;
		File read_file = new File(args.getRef_genome_file());
		BufferedWriter writer = null;
		String line = null;
		StringBuffer qb = new StringBuffer();
		for (int i = 0; i < args.getRead_length(); i ++) {
			qb.append('I');
		}
		String quality = qb.toString();
		try {
			reader = new BufferedReader(new FileReader(read_file));
			writer = new BufferedWriter(new FileWriter(out_file));
			int chr_array = -1;
			boolean first_line = false;
			while((line = reader.readLine()) != null) {
				if (line.charAt(0) == '>') {
					chr_array = Chromosome.chrSymbolToNum(line.substring(1)) - 1;
				}
				else if (chr_array >= 0) {
					ArrayList<Read> chr_read = reads.get(chr_array);
					if (chr_read == null) {
						continue;
					}
					int read_num=0;
					int chr_length = line.length();
					for (int i = 0; i < chr_read.size(); i++) {
						read_num++;
						Read read = chr_read.get(i);
						StringBuffer id = new StringBuffer();
						StringBuffer base_seq = new StringBuffer();
						id.append('@');
						if (read.isCirc_flag()) {
							id.append("circ:");
						}
						if (read.isPeak_flag()) {
							id.append("peak:");
						}
						id.append(read.getId());
						id.append(':');
						id.append(read_num);
						id.append(':');
						id.append(Chromosome.chrNumToSymbol(chr_array + 1));
						id.append(':');
						for (int j = 0; j < read.getPositions().size(); j += 2) {
							id.append(read.getPositions().get(j));
							id.append('-');
							id.append(read.getPositions().get(j + 1));
							id.append('_');
							base_seq.append(line.substring(read.getPositions().get(j) - 1, read.getPositions().get(j + 1)));
						}
						id.setLength(id.length() - 1);
						int base_end = read.getPositions().get(read.getPositions().size() - 1) + 10;
						if (base_end > chr_length) {
							base_end = chr_length;
						}
						base_seq.append(line.substring(read.getPositions().get(read.getPositions().size() - 1), base_end));
						Method.randomErrorSeq(base_seq);
						if (base_seq.length() >= args.getRead_length()) {
							base_seq.setLength(args.getRead_length());
						}
						else {
							System.out.println("Seq too short");
							int to_add = args.getRead_length() - base_seq.length();
							for (int j = 0; j < to_add; j++) {
								base_seq.append('N');
							}
						}
						if (first_line) {
							writer.newLine();
						}
						writer.write(id.toString());
						writer.newLine();
						writer.write(base_seq.toString());
						writer.newLine();
						writer.write('+');
						writer.newLine();
						writer.write(quality);
						first_line = true;
					}
					chr_read.clear();
					writer.flush();
				}
			}
			writer.flush();
			reader.close();
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e1) {
				}
			}
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e1) {
				}
			}
		}
	}
	
	public void readsAppend(ReadInputs args, String out_file, ArrayList<ArrayList<Read>> reads) {
		BufferedReader reader = null;
		File read_file = new File(args.getRef_genome_file());
		BufferedWriter writer = null;
		String line = null;
		StringBuffer qb = new StringBuffer();
		for (int i = 0; i < args.getRead_length(); i ++) {
			qb.append('I');
		}
		String quality = qb.toString();
		try {
			reader = new BufferedReader(new FileReader(read_file));
			writer = new BufferedWriter(new FileWriter(out_file, true));
			int chr_array = -1;
			while((line = reader.readLine()) != null) {
				if (line.charAt(0) == '>') {
					chr_array = Chromosome.chrSymbolToNum(line.substring(1)) - 1;
				}
				else if (chr_array >= 0) {
					ArrayList<Read> chr_read = reads.get(chr_array);
					if (chr_read == null) {
						continue;
					}
					int read_num=0;
					int chr_length = line.length();
					for (int i = 0; i < chr_read.size(); i++) {
						read_num++;
						Read read = chr_read.get(i);
						StringBuffer id = new StringBuffer();
						StringBuffer base_seq = new StringBuffer();
						id.append('@');
						if (read.isCirc_flag()) {
							id.append("circ:");
						}
						if (read.isPeak_flag()) {
							id.append("peak:");
						}
						id.append(read.getId());
						id.append(':');
						id.append(read_num);
						id.append(':');
						id.append(Chromosome.chrNumToSymbol(chr_array + 1));
						id.append(':');
						for (int j = 0; j < read.getPositions().size(); j += 2) {
							id.append(read.getPositions().get(j));
							id.append('-');
							id.append(read.getPositions().get(j + 1));
							id.append('_');
							base_seq.append(line.substring(read.getPositions().get(j) - 1, read.getPositions().get(j + 1)));
						}
						id.setLength(id.length() - 1);
						int base_end = read.getPositions().get(read.getPositions().size() - 1) + 10;
						if (base_end > chr_length) {
							base_end = chr_length;
						}
						base_seq.append(line.substring(read.getPositions().get(read.getPositions().size() - 1), base_end));
						Method.randomErrorSeq(base_seq);
						if (base_seq.length() >= args.getRead_length()) {
							base_seq.setLength(args.getRead_length());
						}
						else {
							System.out.println("Seq too short");
							int to_add = args.getRead_length() - base_seq.length();
							for (int j = 0; j < to_add; j++) {
								base_seq.append('N');
							}
						}
						writer.newLine();
						writer.write(id.toString());
						writer.newLine();
						writer.write(base_seq.toString());
						writer.newLine();
						writer.write('+');
						writer.newLine();
						writer.write(quality);
					}
					chr_read.clear();
					writer.flush();
				}
			}
			writer.flush();
			reader.close();
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e1) {
				}
			}
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e1) {
				}
			}
		}
	}
}
