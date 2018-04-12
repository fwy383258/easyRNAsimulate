package simulator;

import java.io.BufferedWriter;
import java.io.File;
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
	
}
