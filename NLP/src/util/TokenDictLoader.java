package util;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * An object of class TokenDictLoader loads a token list from a file where each
 * line corresponds to a token, and associates a token with its line number
 * (indexing from 0) that is viewed as its token id. This enables one to encode
 * a text into a sequence of numbers according to the tokens' ids defined in the
 * token list file.
 * <p>
 * On the other hand, it is also of interest to translate a sequence of numbers
 * into a text of tokens. We can obtain the mapping between tokens and ids by
 * calling {@link #getIdToTokenMap()}.
 * 
 * @author Youwei Lu
 */
public class TokenDictLoader extends java.util.HashMap<String, Integer> {

	private static final long serialVersionUID = -389245126260027418L;
	private String tokenListFile;
	private String delimiter;
	private int column;

	/**
	 * Loads the map between token strings and token ids from a dictionary file,
	 * with the delimiter set to "\s+" and column set to 0 for
	 * {@link #TokenDictLoader(String, String, int)}
	 * 
	 * @param tokenListFile
	 *            the pathname of the token list file
	 */
	public TokenDictLoader(String tokenListFile) {
		this(tokenListFile, "\\s+", 0);
	}

	/**
	 * Loads the map between token strings and token ids from a dictionary file,
	 * where each row corresponds to a token and the row number (indexing from
	 * zero), strarting from zero, is the token's id.
	 * 
	 * @param tokenListFile
	 *            pathname of the dictionary file
	 * @param delimiter
	 *            a regular expression that use that split rows of the
	 *            dictionary file
	 * @param column
	 *            the column number of each row which stores the token string
	 */
	public TokenDictLoader(String tokenListFile, String delimiter, int column) {
		this.tokenListFile = tokenListFile;
		this.delimiter = delimiter;
		this.column = 0;
	}

	public void loadTokenList() throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(tokenListFile));
		String line, token;
		int i = 0;
		try {
			while ((line = reader.readLine()) != null) {
				token = line.split(delimiter, 2)[column];
				if (put(token, i++) != null) {
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
	 * Returns the map between token ids and token strings, by conversing the
	 * current map between token strings and their ids.
	 * 
	 * @return a map between token ids and token strings
	 */
	public Map<Integer, String> getIdToTokenMap() {
		Map<Integer, String> map = new HashMap<Integer, String>();
		for (Map.Entry<String, Integer> entry : this.entrySet()) {
			map.put(entry.getValue(), entry.getKey());
		}
		return map;
	}
}
