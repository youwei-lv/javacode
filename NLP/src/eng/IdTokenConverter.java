package eng;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.SProperties;
import util.TokenDictLoader;

/**
 * An object of class IdTokenConvert translates token ids contained in files to
 * their original textual forms defined in a token list. An IdTokenConverter
 * scans each row of an input file, extract the string containing token ids
 * according to the regular expression characterizing the format of each line
 * and besides the group number of patterns for the token ids, and finally
 * interpolate the line by replacing the token ids with the token strings to
 * obtain the translated line of text. The results will be written to an output
 * file.
 * <p>
 * 
 * There are several arguments that need passing to the command line.
 * <ul>
 * <li>tokenlist.file - the pathname of the token list file</li>
 * <li>tokenlist.delimiter - the delimiter used in the token list file to
 * separate the content of each row</li>
 * <li>tokenlist.column - the column number of the token string in each row,
 * with indexing from zero</li>
 * <li>lineStringPattern - a regular expression for characterizing the format of
 * each line in token id files, which may consist of several groups, and the
 * token ids can be identified by a given group number</li>
 * <li>groupOfTokenId - the group number of token id string in the
 * lineStringPattern, with indexing from 1.</li>
 * <li>input - the input which can be a normal file or a directory</li>
 * <li>output - the output which stands for the output file provided that the
 * input file is a normal file; otherwise, it means a directory</li>
 * <li>showWhitespace - indicates whether or not the whitespace between token
 * strings would be written to the output file</li>
 * </ul>
 * 
 * @author Youwei Lu
 */
public class IdTokenConverter {
	private String lineStringPattern;
	private int groupOfTokenId;
	private TokenDictLoader dictLoader;
	private Map<Integer, String> idTokenMap;
	private String inputf;
	private String outputf;
	private boolean isDir;
	private boolean dispWhitespace;
	private Pattern linePattern;

	public IdTokenConverter(String tokenListFile, String delimiter, int column, String lineStringPattern,
			int groupOfTokenId, String inputf, String outputf, boolean dispWhitespace) throws IOException {
		File inputfile = new File(inputf);
		File outputfile = new File(outputf);
		if (!inputfile.exists()) {
			throw new RuntimeException(inputf + " does not exist.");
		}
		if (inputfile.isDirectory()) {
			if (outputfile.exists() && outputfile.isFile()) {
				throw new RuntimeException(outputf + " cannot be an existing normal file.");
			} else if (!outputfile.exists()) {
				outputfile.mkdirs();
			}
			this.isDir = true;
		} else {
			this.isDir = false;
		}
		this.dictLoader = new TokenDictLoader(tokenListFile, delimiter, column);
		this.lineStringPattern = lineStringPattern;
		this.groupOfTokenId = groupOfTokenId;
		this.inputf = inputf;
		this.outputf = outputf;
		this.dispWhitespace = dispWhitespace;
		dictLoader.loadTokenList();
		idTokenMap = dictLoader.getIdToTokenMap();
		linePattern = Pattern.compile(this.lineStringPattern);
	}

	private void translateTokenFile(File idFile, File tokenFile) throws IOException {
		BufferedReader read = new BufferedReader(new FileReader(idFile));
		BufferedWriter write = new BufferedWriter(new FileWriter(tokenFile));
		String line = null, leftPart, rightPart, idPart, idString;
		int p1, p2, id;
		Matcher m;
		Pattern idPattern = Pattern.compile("\\d+");
		StringBuilder buffer = null;

		while ((line = read.readLine()) != null) {
			/* extract token-id and non-token-id strings */
			m = linePattern.matcher(line);
			if (!m.find()) {
				continue;
			}
			buffer = new StringBuilder();
			idPart = m.group(groupOfTokenId);
			leftPart = line.substring(0, m.start(groupOfTokenId));
			rightPart = line.substring(m.end(groupOfTokenId));
			m = idPattern.matcher(idPart);
			p1 = 0;
			while (m.find()) {
				p2 = m.start();
				if (!idPart.substring(p1, p2).matches("\\s+") || this.dispWhitespace) {
					buffer.append(idPart.substring(p1, p2));
				}
				p1 = m.end();
				idString = idPart.substring(p2, p1);
				id = Integer.valueOf(idString);
				if (idTokenMap.containsKey(id)) {
					buffer.append(idTokenMap.get(id));
				} else {
					buffer.append(idString);
				}
			}
			buffer.append(idPart.substring(p1));
			write.write(leftPart + buffer.toString() + rightPart + "\n");
		}
		write.flush();
		write.close();
		read.close();
	}

	public void doTranslate() throws IOException {
		File input = new File(inputf);
		File output = new File(outputf);

		if (isDir) {
			for (File f : input.listFiles()) {
				translateTokenFile(new File(input, f.getName()), new File(output, f.getName()));
			}
		} else {
			translateTokenFile(input, output);
		}
	}

	public static void main(String[] args) throws IOException {
		SProperties props = new SProperties();
		props.load(args);

		String tokenListFile = "/home/youwei/Github/npbayes-hsmm/unif-thsmm-unit_pure_C/dict.txt";
		String delimiter = "\\s+";
		int column = 0;
		String lineStringPattern = "(.+):(.+)";
		int groupOfTokenId = 2;
		String inputf = "/home/youwei/Github/npbayes-hsmm/unif-thsmm-unit_pure_C/wt/0.1";
		String outputf = "/home/youwei/iiiiiii";
		boolean dispWhitespace = false;

		tokenListFile = props.getProperty("tokenlist.file");
		delimiter = props.getProperty("tokenlist.delimiter");
		column = Integer.valueOf(props.getProperty("tokenlist.column"));
		lineStringPattern = props.getProperty("lineStringPattern");
		groupOfTokenId = Integer.valueOf(props.getProperty("groupOfTokenId"));
		inputf = props.getProperty("input");
		outputf = props.getProperty("output");
		dispWhitespace = Boolean.valueOf(props.getProperty("showWhitespace"));

		IdTokenConverter converter = new IdTokenConverter(tokenListFile, delimiter, column, lineStringPattern,
				groupOfTokenId, inputf, outputf, dispWhitespace);
		converter.doTranslate();
	}
}
