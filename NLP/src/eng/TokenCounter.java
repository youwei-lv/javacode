package eng;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

/**
 * A TokenCounter object counts the appearance times of different tokens in
 * files; however, there are two regular expressions for token selection, one
 * for collecting tokens and the other for scraping tokens, and the counted
 * tokens are those that conform to the {@linkplain #countPattern collecting
 * regular expression} and that do not match the {@linkplain #ignorePattern
 * scraping one}.
 * 
 * @author Youwei Lu
 */
public class TokenCounter {

	/**
	 * Specified counting target, which may be a file name or a directory. When
	 * it is a single file, the tokens in the file will be counted; however,
	 * when it is a directory, tokens in all files in the directory will be
	 * counted.
	 */
	String inputf;
	/**
	 * The name of file to which the counting result will be written.
	 */
	String countf;
	/**
	 * A regular expression for matched tokens to be filtered out from counting
	 */
	String ignorePattern;
	/**
	 * A regular expression for matched tokens to be counted
	 */
	String countPattern;
	/**
	 * A counting map between tokens and their frequencies
	 */
	private Map<String, Integer> count;
	/**
	 * Constant standing for sorting in ascending order.
	 */
	static public final int ASC = 0;
	/**
	 * Constant standing for sorting in descending order.
	 */
	static public final int DESC = 1;

	public TokenCounter(String inputf, String countf) {
		this(inputf, countf, ".+", "");
	}

	private TokenCounter(String inputf, String countf, String countPattern, String ignorePattern) {
		this.inputf = inputf;
		this.countf = countf;
		this.countPattern = countPattern;
		this.ignorePattern = ignorePattern;
		count = new HashMap<String, Integer>();
	}

	/**
	 * Increase the frequency of a specific token by one, if the token matches
	 * the collecting regular expression and does not match the scraping regular
	 * expression.
	 */
	public void addToken(String token) {
		if (Pattern.matches(this.ignorePattern, token) == false && Pattern.matches(this.countPattern, token) == true) {
			if (count.containsKey(token)) {
				count.put(token, count.get(token) + 1);
			} else {
				count.put(token, 1);
			}
		}
	}

	/**
	 * Output a list of token-frequency pairs which are sorted in ascending or
	 * descending order, based on tokens' frequency.
	 * 
	 * @param order
	 *            Mode of sorting, which is either {@link #ASC} or {@link #DESC}
	 * @return A sorted list of token-frequency pairs.
	 */
	public List<Map.Entry<String, Integer>> sortByFrequency(int order) {
		List<Map.Entry<String, Integer>> list = new LinkedList<>(count.entrySet());
		if (order == ASC) {
			Collections.sort(list, new Comparator<Map.Entry<String, Integer>>() {
				@Override
				public int compare(Map.Entry<String, Integer> o1, Map.Entry<String, Integer> o2) {
					return o1.getValue().compareTo(o2.getValue());
				}
			});
		} else if (order == DESC) {
			Collections.sort(list, new Comparator<Map.Entry<String, Integer>>() {
				@Override
				public int compare(Map.Entry<String, Integer> o1, Map.Entry<String, Integer> o2) {
					return -o1.getValue().compareTo(o2.getValue());
				}
			});
		} else {
			throw new RuntimeException("Unknown sorting order.");
		}
		return list;
	}

	public void docount() {
		File input = new File(this.inputf);
		BufferedReader reader = null;
		String line = null;
		if (input.isDirectory()) {
			File[] directoryListing = input.listFiles();
			for (File child : directoryListing) {
				System.out.println("processing file " + child.getName());
				try {
					reader = new BufferedReader(new FileReader(child));
					while ((line = reader.readLine()) != null) {
						String tokens[] = line.split("\\s+");
						for (String token : tokens) {
							addToken(token);
						}
					}
				} catch (IOException e) {
					e.printStackTrace();
				} finally {
					if (reader != null)
						try {
							reader.close();
						} catch (IOException e) {
							e.printStackTrace();
						}
				}
			}
		} else {
			try {
				reader = new BufferedReader(new FileReader(input));
				while ((line = reader.readLine()) != null) {
					String tokens[] = line.split("\\s+");
					for (String token : tokens) {
						addToken(token);
					}
				}
			} catch (IOException e) {
				e.printStackTrace();
			} finally {
				if (reader != null)
					try {
						reader.close();
					} catch (IOException e) {
						e.printStackTrace();
					}
			}
		}
	}

	public void outputCountMap(List<Map.Entry<String, Integer>> list) throws IOException {
		BufferedWriter writer = new BufferedWriter(new FileWriter(this.countf));
		for (Map.Entry<String, Integer> entry : list) {
			writer.write(entry.getKey() + " " + entry.getValue());
			writer.newLine();
		}
		writer.close();
	}

	public String getIgnorePattern() {
		return ignorePattern;
	}

	public void setIgnorePattern(String ignorePattern) {
		this.ignorePattern = ignorePattern;
	}

	public String getCountPattern() {
		return countPattern;
	}

	public void setCountPattern(String countPattern) {
		this.countPattern = countPattern;
	}

	public static void main(String[] args) throws IOException {
		String input, output, countPattern, ignorePattern;
		TokenCounter counter;
		if ((input = System.getProperty("input")) == null) {
			throw new RuntimeException("Use -Dinput=<file name> option to specify the input file/directory.");
		}
		if ((output = System.getProperty("output")) == null) {
			throw new RuntimeException("Use -Doutput=<file name> option to specify the output file.");
		}

		counter = new TokenCounter(input, output);

		if ((countPattern = System.getProperty("countpattern")) != null) {
			counter.setCountPattern(countPattern);
		}

		if ((ignorePattern = System.getProperty("ignorepattern")) != null) {
			counter.setIgnorePattern(ignorePattern);
		}

		counter.docount();
		counter.outputCountMap(counter.sortByFrequency(TokenCounter.DESC));
	}

}