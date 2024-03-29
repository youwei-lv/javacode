package eng;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Pattern;

import util.SProperties;

/**
 * A TokenIdConverter convert the tokens in a single file or the files in a
 * directory to their ids which are separated by a delimiter "|" which stands for
 * the separators in the original text, e.g. commas, newline characters etc., to
 * identify segments in the resulted id sequence; the id of a token is its line
 * number (indexing from zero) in the token list file where each line begins with
 * a token.
 * <p>
 * 
 * Because some words in the files can be missing from the token list,
 * extracting all separators in a file may result in some empty segments. This
 * situation is considered and is handled by automatically removing these empty
 * segments.
 * <p>
 * 
 * When the input is a single text file, the output will also be a text file. If
 * the specified output is an existing directory, then an error will occur;
 * however, when the input is a directory consisting of multiple text files, the
 * specified output stands for a directory to take in the resulted files. For
 * each token file in the input directory, a file with the same filename will be
 * created in the output directory to represent the resulted id sequence by
 * processing the token file.
 * <p>
 * 
 * Compiled TokenIdConverter can be called through command line as a tool to
 * convert files of tokens into files of ids according to a token list file. To
 * specify these related files, use the following options:
 * <ul>
 * <li>-input: denotes the input of files (a single file or a directory of
 * files)</li>
 * <li>-output: denotes the output of files (a single file or a directory of
 * files)</li>
 * <li>-tokenlist: denotes the filename of the token list file, where each line
 * begins with token string, and the id is the line number of the token</li>
 * <li>-printword: a switch determining whether the token are shown following
 * the converted ids</li>
 * <li>-segmentpattern: the regular expression for matching the segments in the
 * origin texts, by default set to \\.|,|\\n|:|\\?|!</li>
 * </ul>
 * An example of specified options: -input=test_output -output=test_output_ids
 * -printword=false -tokenlist=token_freq -printword=true
 * 
 * @author Youwei Lu
 */
public class TokenIdConverter {

	/**
	 * The filename of the token list file
	 */
	private String tokenListFile;
	/**
	 * The map between the token and its id. A token's id is the line number
	 * whose corresponding line begins with the token
	 */
	private Map<String, Integer> tokenIdMap;
	/**
	 * A regular expression for matching the separators in files
	 */
	private String separatorPattern;

	public TokenIdConverter(String tokenListFile) {
		this.tokenListFile = tokenListFile;
		this.separatorPattern = "\\.|,|\\n|:|\\?|!";
		tokenIdMap = new HashMap<String, Integer>();
	}

	public void loadTokenList() throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(tokenListFile));
		String line, token;
		int i = 0;
		try {
			while ((line = reader.readLine()) != null) {
				token = line.split("\\s+", 2)[0];
				if (tokenIdMap.put(token, i++) != null) {
					throw new RuntimeException(
							"Token " + token + " has multiple ids in token list file " + tokenListFile);
				}
			}
		} catch (IOException e) {
			throw e;
		} finally {
			if (reader != null)
				reader.close();
		}
	}

	/**
	 * Read tokens from an input file, convert the tokens to corresponding ids
	 * in the token list, and write the results to the output file. If
	 * separators are detected, then insert a delimiter "|" in the id sequence
	 * to denote segments; however, empty segments are avoided.
	 * 
	 * @param infile
	 *            a File object representing the input file
	 * @param outfile
	 *            a File object representing the output file
	 * @param printword
	 *            an indicator for printing the token or not
	 * @throws java.io.IOException
	 *             thrown when infile does not exist or the writing to outfile
	 *             fails
	 */
	private void convertTokenFile(File infile, File outfile, boolean printword) throws IOException {
		String line, tokens[], token;
		int prevId = -1, id = -1, i, j;
		BufferedReader read = null;
		BufferedWriter write = null;

		try {
			read = new BufferedReader(new FileReader(infile));
			write = new BufferedWriter(new FileWriter(outfile));

			while ((line = read.readLine()) != null) {
				// readLine() removes the ending newline character, so append it manually
				line = line + "\n";
				tokens = line.split("( |\\t|\\r)+");
				for (i = 0; i < tokens.length; i++) {
					token = tokens[i];
					j = i + 1;
					if (Pattern.matches(separatorPattern, token)) {
						/*
						 * if the current token is a separator, it is written to
						 * the output only when it is between two ids
						 */
						for (; j < tokens.length; j++) {
							if (tokenIdMap.containsKey(tokens[j])) {
								prevId = id;
								id = tokenIdMap.get(tokens[j]);
								break;
							}
						}
						if (j < tokens.length) {
							if (prevId >= 0) {
								if (printword) {
									write.write("| " + id + "(" + tokens[j] + ")" + " ");
								} else {
									write.write("| " + id + " ");
								}
							} else {
								if (printword) {
									write.write(id + "(" + tokens[j] + ")" + " ");
								} else {
									write.write(id + " ");
								}
							}
						} else if (prevId >= 0 && j == tokens.length) {
							write.write("| ");
						}
						i = j;
					} else if (tokenIdMap.containsKey(token)) {
						prevId = id;
						id = tokenIdMap.get(token);
						if (printword) {
							write.write(id + "(" + token + ")" + " ");
						} else {
							write.write(id + " ");
						}
					} else {
						prevId = id;
						id = -1;
					}
				}
			}
		} catch (IOException e) {
			throw e;
		} finally {
			write.close();
			read.close();
		}
	}

	public String getSeparatorPattern() {
		return separatorPattern;
	}

	public void setSeparatorPattern(String separatorPattern) {
		this.separatorPattern = separatorPattern;
	}

	private void convertTokenDirectory(File inputdir, File outputdir, boolean printword) throws IOException {
		File[] directoryListing = inputdir.listFiles();
		File inputfile, outputfile;

		for (File child : directoryListing) {
			System.out.println("processing file " + child.getName());
			inputfile = new File(inputdir, child.getName());
			outputfile = new File(outputdir, child.getName());
			convertTokenFile(inputfile, outputfile, printword);
		}
	}

	public void doconvert(File input, File output, boolean printword) throws IOException {
		if (input.isFile()) {
			if (output.isDirectory()) {
				throw new RuntimeException(
						"Error: " + input + "is a single file, while" + output + " is an existing directory.");
			}
			convertTokenFile(input, output, printword);
		} else {
			if (output.isFile()) {
				throw new RuntimeException(
						"Error: " + input + " is a directory, while " + output + " is an existing file.");
			}
			if (!output.exists()) {
				output.mkdirs();
			}
			convertTokenDirectory(input, output, printword);
		}
	}

	public static void main(String[] args) throws IOException {
		SProperties props = new SProperties();
		props.load(args);
		File inputf = new File(props.getProperty("input"));
		File outputf = new File(props.getProperty("output"));
		String tokenList = props.getProperty("tokenlist");
		boolean printword = Boolean.parseBoolean(props.getProperty("printword"));

		TokenIdConverter convert = new TokenIdConverter(tokenList);
		if (props.containsKey("segmentpattern")) {
			convert.setSeparatorPattern(props.getProperty("segmentpattern"));
		}
		convert.loadTokenList();
		convert.doconvert(inputf, outputf, printword);
	}

}
