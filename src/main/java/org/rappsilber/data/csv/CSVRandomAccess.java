/*
 * Copyright 2015 Lutz Fischer <lfischer at staffmail.ed.ac.uk>.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.rappsilber.data.csv;

import java.awt.event.ActionEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.URI;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.rappsilber.utils.UpdatableChar;

/**
 * A CSV-parser, that will cache the csv-file in memory and permits a access to the data by row and column
 * @author lfischer
 */
public class CSVRandomAccess extends CsvParser {
    ArrayList<String[]> m_data = new ArrayList<String[]>();;
    int m_current = -1;
    Boolean m_loading = false;
    private final static Pattern NONNUMERIC = Pattern.compile("[^0-9E\\-.+]*");    
    private final static Pattern NUMERIC = Pattern.compile("[^0-9]*([\\+\\-]?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:\\s?E[\\+\\s]?[0-9]+)?).*");    
    
    /**
     * CSVCondition is used to find rows in a csv-file that have a specific value in the given column
     */
    public class CSVCondition {
        int field;
        String value;
        
        /**
         * creates a new condition
         * @param field
         * @param value
         */
        public CSVCondition(int field, String value) {
            this.field = field;
            this.value = value;
        }
        
        /**
         * returns whether the given row fits the condition
         * @param row
         * @return true if it does; false otherwise
         */
        boolean fits(int row) {
            return m_data.get(row)[field].contentEquals(value);
        } 
        
    }
    
    /**
     * Enum with constants for defining a how to sort a file
     */
    public static enum CSVSortType {

        /**
         * sort numeric - fileds are treated as double values and sorted by value
         */
        numeric,

        /**
         * sorted as string values
         */
        alphanumeric
    }
    
    /**
     * base class that is used for comparing rows during sorting
     */
    public static abstract class CSVSort {
        int field;

        /**
         *
         * @param field
         */
        public CSVSort(int field) {
            this.field = field;
        }
        
        abstract int compare(String[] row1, String[] row2);
        
    }

    /**
     * implement a numeric compare of a given field for sorting the file
     */
    public static class CSVSortNumeric extends CSVSort {

        /**
         *
         * @param field
         */
        public CSVSortNumeric(int field) {
            super(field);
        }
        
        int compare(String[] row1, String[] row2) {
            Double d1;
            Double d2;
            try {
                d1 = Double.parseDouble(row1[field]);
            } catch (Exception nfe) {
                d1 = Double.NEGATIVE_INFINITY;
            } 
            try {
                d2 = Double.parseDouble(row2[field]);
            } catch (Exception nfe) {
                d2 = Double.NEGATIVE_INFINITY;
            }
            return d1.compareTo(d2);
        }
        
    }

    /**
     * implement a string compare of a given field for sorting the file
     */
    public static class CSVSortAlphaNumeric extends CSVSort {

        /**
         * constructor
         * @param field by what column to sort
         */
        public CSVSortAlphaNumeric(int field) {
            super(field);
        }
        
        @Override
        int compare(String[] row1, String[] row2) {
            return row1[field].compareTo(row2[field]);
        }
        
    }
    
    /**
     * The loading and processing can be done asynchronous - listeners can be used to get information on completions 
     */
    public static interface LoadListener {

        /**
         * function to be called when the assigned event fires
         * @param row
         * @param column
         */
        void listen(int row,int column);
    };
    
    /**
     * list of listeners to be notified if a file loading was finished
     */
    ArrayList<LoadListener> m_listenerCompleted = new ArrayList<LoadListener>();
    /**
     * list of listeners to be notified regularly for the progress of loading the file
     */
    ArrayList<LoadListener> m_listenerProgress = new ArrayList<LoadListener>();
    /**
     * Listeners to be notified when sorting is finished
     */
    ArrayList<LoadListener> m_listenerSort = new ArrayList<LoadListener>();
    /**
     * Listeners to be notified when a colum was added, removed or renamed
     */
    ArrayList<LoadListener> m_listenerColumnsChanged = new ArrayList<LoadListener>();
    /**
     * Listeners to be notified when the header of the csv-file was modified
     */
    ArrayList<java.awt.event.ActionListener> m_listenerHeaderChanged = new ArrayList<java.awt.event.ActionListener>();
    

    /**
     * Creates a new instance - without reading any file. 
     * @param delimiter the delimiter to be used
     * @param quote  the quote character to be used
     */
    public CSVRandomAccess(char delimiter, char quote) {
        super(delimiter, quote);
    }

    /**
     * Creates a new instance - and interpretes the provides String as CSV data to be read in
     * @param delimiter
     * @param quote 
     * @param data the data to interpreted as csv-file
     */
    public CSVRandomAccess(char delimiter, char quote, String data) {
        this(delimiter, quote);
        String[] s = data.split("\n");
        int mc = 0;
        for (String line : s) {
            ArrayList<String> sl = splitLine(line);
            String[] sla = new String[sl.size()];
            sla = sl.toArray(sla);
            m_data.add(sla);
            if (sla.length > mc)
                mc = sla.length;
        }
        setMaxColumns(mc);
    }
    
    /**
     * Creates a new instance - without reading any file. 
     */
    private CSVRandomAccess() {
        super();
    }

    /**
     * returns a condition, that can be used to search for rows that have the given value
     * @param field column where to look for the value
     * @param value value to look for
     * @return a new CSVConditin
     */
    public CSVCondition getCondition(int field, String value) {
        return new CSVCondition(field, value);
    }
    
    /**
     * Read in a file as CSV-file
     * @param f the file to be read
     * @param hasHeader is the first row in the file to be treated as header (true) or first data row (false)
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public void openFile(File f, boolean hasHeader) throws FileNotFoundException, IOException {
        super.openFile(f,false);
        m_loading = true;
        int row = 0;
        synchronized(m_data) {
            if (hasHeader && super.next()) {
                m_current = 0;
                m_data.add(super.getValues());
                setCurrentLineAsHeader();
                notifyProgress(row);
            }
            while (super.next()) {
                if (row++ % 10 == 0)
                    notifyProgress(row);
                m_data.add(super.getValues());
            }
            m_loading = false;
        }
//        m_loading.notifyAll();
        notifyComplete();
    }

    
    public void openURI(URI f, boolean hasHeader) throws FileNotFoundException, IOException {
        if (f.getScheme().equals("file")) {
            openFile(new File(f), hasHeader);
        } else {
            super.openURI(f, false);
            m_loading = true;
            int row = 0;
            synchronized(m_data) {
                if (hasHeader && super.next()) {
                    m_current = 0;
                    m_data.add(super.getValues());
                    setCurrentLineAsHeader();
                    notifyProgress(row);
                }
                while (super.next()) {
                    if (row++ % 10 == 0)
                        notifyProgress(row);
                    m_data.add(super.getValues());
                }
                m_loading = false;
            }
    //        m_loading.notifyAll();
            notifyComplete();            
        }
    }
    
    public void openURI(URI f) throws FileNotFoundException, IOException {
        openURI(f, false);
    }
    
    /**
     * Read in a file as CSV-file - this function returns as soon as the file is open - and if  hasHeader is true - the header-row was read in.
     * <br/> To be notified when the file is completely read in the {@link #addListenerComplete(org.rappsilber.data.csv.CSVRandomAccess.LoadListener) addListenerComplete} can be used.
     * <br/> Intermediate Notifications can be received via {@link #addListenerProgress(org.rappsilber.data.csv.CSVRandomAccess.LoadListener)  addListenerProgress}.
     * @param f the file to be read
     * @param hasHeader is the first row in the file to be treated as header (true) or first data row (false)
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public void openFileAsync(File f, boolean hasHeader) throws FileNotFoundException, IOException {
        super.openFile(f,false);
        m_loading = true;
        
        // read the first line as header if requested and a line is there
        if (hasHeader && super.next()) {
            m_current = 0;
            m_data.add(super.getValues());
            setCurrentLineAsHeader();
            notifyProgress(0);
        }
        
        Runnable runnable = new Runnable() {

            @Override
            public void run() {
                try {
                    synchronized(m_data) {
                        int row = 0;
                        while (CSVRandomAccess.super.next()) {
                            if (row++ % 10 == 0) {
                                notifyProgress(row);
                            }
                            m_data.add(CSVRandomAccess.super.getValues());
                        }
                        m_loading = false;
                    }
        //        m_loading.notifyAll();
                    notifyComplete();

                } catch (IOException ex) {
                    Logger.getLogger(CSVRandomAccess.class.getName()).log(Level.SEVERE, null, ex);
                    notifyComplete();
                }
            }
        };
        new Thread(runnable).start();
    }
    
    /**
     * move the active-row-pointer to the next row.
     * @return true if there is actually a next row - or false if there is no next row
     */
    @Override
    public boolean next() {
        if (m_current<m_data.size() -1) {
            m_current ++;
            setCurrentValues(m_data.get(m_current));
            return true;
        }
        m_current = m_data.size();
        return false;
    }
    
    /**
     * move the active-row-pointer to the previous row.
     * @return true if the the last active row is not the first row.
     */
    public boolean previous() {
        if (m_current>0) {
            m_current --;
            setCurrentValues(m_data.get(m_current));
            return true;
        }
        m_current = -1;
        return false;
    }

    /**
     * is the "current" row within the csv-file
     * @return 
     */
    public boolean activeRow() {
        return m_current >=0 && m_current< m_data.size();
    }
    
    /**
     * set the given row as the current row
     * @param row 
     */
    public void setRow(int row) {
        m_current = row -1;
        next();
    }
    
    /**
     * removes the currently active row
     */
    public void deleteCurrentLine() {
        synchronized(m_data) {
            deleteRow(m_current); 
        }
        
    }
    
    /**
     * removes the specified row
     * @param row row to be deleted
     */
    public void deleteRow(int row) {
        synchronized(m_data) {
            m_data.remove(row);
            if (row >= m_current)
                m_current--;
        }
    }
    
    /**
     * find the first row, that matches the condition. If a row is found that row is then defined as the current row
     * @param conditions
     * @return the first row that matches, if no row matches then -1 is returned
     */
    public int findFirstRow(CSVCondition[] conditions) {
        row: for (int i = 0; i<m_data.size(); i++) {
            for (int c =0; c< conditions.length; c++) {
                if (!conditions[c].fits(i))
                    continue row;
            }
            m_current = i;
            return i;
        }
        return -1;
    }

    /**
     * find the next row after the current row, that matches the condition. 
     * If a row is found that row is then defined as the current row
     * @param conditions
     * @return the first row that matches, if no row matches then -1 is returned
     */
    public int findFromCurrent(CSVCondition[] conditions) {
        
        row: for (int i = Math.max(m_current+1,0); i<m_data.size(); i++) {
            for (int c =0; c< conditions.length; c++) {
                if (!conditions[c].fits(i))
                    continue row;
            }
            m_current = i;
            return i;
        }
        return -1;
    }

    /**
     * creates a new CSVRandomAccess instance, that has the same data as the current instance.
     * The data are semi cloned. each instance has a its own list of rows - but the data within these rows are shared.
     * In effect if a row is deleted in on Instance it is not deleted in the second. But if a value is changed in one the value will also be changed in the second instance.
     * @return the cloned instance
     */
    @Override
    public CSVRandomAccess  clone() {
        CSVRandomAccess ccsv = new CSVRandomAccess(getDelimiter().charAt(0), getQuote().charAt(0));
        super.transfer(ccsv);
        synchronized(m_data) {
            ccsv.m_data  = (ArrayList<String[]>) m_data.clone();
        }
        ccsv.m_current  = m_current;
        ccsv.m_loading  = m_loading;
        return ccsv;
    }
    
    /**
     * takes the current line, and sets it as header
     */    
    public void setCurrentLineAsHeader() {
        setHeader(getCurrentLine());
        deleteCurrentLine();
        fireHeaderChanged();
    }
    
    /**
     * takes the given line, and sets it as header
     * @param row - row to be used as header
     */    
    public void setLineAsHeader(int row) {
        int oldCurrentlin = m_current;
        m_current = row;
        setCurrentLineAsHeader();
        m_current = oldCurrentlin;
//        fireHeaderChanged();
    }

    /**
     * inserts the given line as a new row in the csv-file (in memory only)
     * @param row where to insert the row
     * @param line the row to be inserted
     */    
    public void insertLine(int row,String line) {        
        ArrayList<String> sl = splitLine(line);
        String[] d = new String[sl.size()];
        d = sl.toArray(d);
        insertLine(row, d);
    }

    /**
     * inserts the given line as a new row in the csv-file (in memory only)
     * @param row where to insert the row
     * @param line the row to be inserted
     */    
    public void insertLine(int row,String[] line) {        
        m_data.add(row, line);
        if (line.length > getMaxColumns())
            setMaxColumns(line.length);
        if (m_current >= row)
            m_current++;
    }
    
    /**
     * Creates a new instance of CSVRandomAccess for a given file while trying 
     * to guess the right values for the delimeter and quote characters
     * @param f the file to be read in
     * @param guessLines - how many lines to read in to guess the delimiter and quote characters
     * @param hasHeader - does the file have a header - row
     * @param delimiters list of delimiters to test
     * @param quotes list of quote characters to test
     * @return a new instance of CSVRandomAccess
     * @throws FileNotFoundException 
     * @throws IOException 
     */
    public static CSVRandomAccess guessCsvRA(File f, int guessLines, boolean hasHeader, char[] delimiters, char[] quotes) throws FileNotFoundException, IOException {
      
        UpdatableChar delimiter = new UpdatableChar();
        UpdatableChar quote = new UpdatableChar();
        Boolean unique = guessDelimQuote(f, guessLines, delimiters, quotes, delimiter, quote);
                
        
        if (unique == null) {
            return null;
        }
        
        CSVRandomAccess csv = new CSVRandomAccess(delimiter.value, quote.value);
        csv.openFile(f, hasHeader);
        return csv;
    }

    /**
     * Creates a new instance of CSVRandomAccess for a given file while trying 
     * to guess the right values for the delimeter and quote characters.
     * The actual reading will be done in a separate thread and the function will return immediately
     * @param f the file to be read in
     * @param guessLines - how many lines to read in to guess the delimiter and quote characters
     * @param hasHeader - does the file have a header - row
     * @param delimiters list of delimiters to test
     * @param quotes list of quote characters to test
     * @return a new instance of CSVRandomAccess
     * @throws FileNotFoundException 
     * @throws IOException 
     */
    public static CSVRandomAccess guessCsvAsync(final File f, int guessLines, final boolean hasHeader, char[] delimiters, char[] quotes) throws FileNotFoundException, IOException {
      
        UpdatableChar delimiter = new UpdatableChar();
        UpdatableChar quote = new UpdatableChar();
        Boolean unique = guessDelimQuote(f, guessLines, delimiters, quotes, delimiter, quote);
                
        
        if (unique == null) {
            return null;
        }
        
        final CSVRandomAccess csv = new CSVRandomAccess(delimiter.value, quote.value);
        csv.openFileAsync(f, hasHeader);

        return csv;
    }
    
    /**
     * Opens a csv-parser for the given file and tries to predict the quote 
     * character and the field-delimiter based on the first number of lines.<br/>
     * After it reads the first lines it will reopen the file with a configured 
     * csv-parser.<br/>
     * <b><i>This means, that the file must be re-readable e.g no pipe</i></b>
     * 
     * @param f the file to open
     * @param guessLines how many lines to read in for guessing the parameters
     * @param hasHeader - does the file have a header - row
     * @return a {@link CSVRandomAccess} for the file
     * @throws FileNotFoundException
     * @throws IOException
     */
    public static CSVRandomAccess guessCsvRA(File f, int guessLines, boolean hasHeader) throws FileNotFoundException, IOException {
        return guessCsvRA(f, guessLines, hasHeader, TEST_DELIMITERS, TEST_QUOTES);
    }

    /**
     * Opens a csv-parser for the given file and tries to predict the quote 
     * character and the field-delimiter based on the first number of lines.<br/>
     * After it reads the first lines it will reopen the file with a configured 
     * csv-parser.<br/>
     * <b><i>This means, that the file must be re-readable e.g. no pipe.</i></b>
     * <p>
     * This will read the file asynchronously. Meaning the function will return after the initial guessing. 
     * </p>
     * @param f the file to open
     * @param guessLines how many lines to read in for guessing the parameters
     * @param hasHeader - does the file have a header - row
     * @return a {@link CSVRandomAccess} for the file
     * @throws FileNotFoundException
     * @throws IOException
     */
    public static CSVRandomAccess guessCsvAsync(File f, int guessLines, boolean hasHeader) throws FileNotFoundException, IOException {
        return guessCsvAsync(f, guessLines, hasHeader, TEST_DELIMITERS, TEST_QUOTES);
    }
    
    /**
     * returns the specified value. 
     * If the field does not contain a value or defined field is outside the CSV-file {@link #MISSING_FIELD} is returned
     * @param field column 
     * @param row row
     * @return value at the specified place.
     */
    public String getValue(String field, Integer row) {
        if (row ==null || row >= m_data.size())
            return MISSING_FIELD;
        
        Integer column = getColumn(field);
        
        if (column == null)
            return MISSING_FIELD;
        
        String[] line = m_data.get(row); 
        if (column >= line.length)
            return MISSING_FIELD;
        
        return line[column];
    }
    
    
    /**
     * returns the specified value. 
     * If the field does not contain a value or defined field is outside the CSV-file defaultValue is used to calculate a return value
     * @param field column 
     * @param row row
     * @param defaultValue provides a default value for each row
     * @return value at the specified place.
     */
    public String getValue(Integer field, Integer row, CSVValueCalc defaultValue) {
        if (field == null || row == null || row >= m_data.size())
            return defaultValue.getValue(this, row);

        String[] line = m_data.get(row); 

        if (field >= line.length)
            return defaultValue.getValue(this, row);

        // we have an entry for it but it is MISSING_FIELD -> so we should just replace it
        if (line[field] == MISSING_FIELD) {
            String v = defaultValue.getValue(this);
            setValue(v, field, row);
            return v;
        }
            
        return line[field];
    }
    
    
    /**
     * returns the specified value. 
     * If the field does not contain a value or defined field is outside the CSV-file defaultValue is used to calculate a return value
     * @param field column 
     * @param row row
     * @param defaultValue provides a default value for each row
     * @return value at the specified place.
     */
    public String getValue(String field, Integer row, CSVValueCalc defaultValue) {
        Integer column = getColumn(field);
        if (column == null)
            return defaultValue.getValue(this,row);
        return getValue(column, row, defaultValue);
    }        
    
    /**
     * Sets the value of the given field as value.toString()
     * <br/>
     * If the row did not previously had enough columns to accommodate the field,
     * the needed amount of fields will be added with MISSING_VALUE.
     * 
     * @param value the new value
     * @param field column of the 
     * @param row   row where the value should be changed
     */
    public void setValue(Object value, int field, int row) {
        if (row >= m_data.size()) {
            if (!value.equals(MISSING_FIELD)) {
                synchronized (m_data) {
                    for (int i = m_data.size(); i < row; i++) {
                        m_data.add(new String[0]);
                    }
                }
                String[] line = new String[field+1];
                for (int i = 0; i<field; i++) {
                    line[i] = MISSING_FIELD;
                }
                line[field] = value.toString();
                m_data.add(line);
                if (field >= getMaxColumns())
                    setMaxColumns(field + 1);
            }
        } else {
            String[] line = m_data.get(row); 
            if (field >= line.length) {
                if (!value.equals(MISSING_FIELD)) {
                    String[] dummyField = java.util.Arrays.copyOf(line, field + 1);
                    for (int f = line.length; f< field; f++) {
                        dummyField[f] = MISSING_FIELD;
                    }
                    dummyField[field] = value.toString();
                    m_data.set(row, dummyField);
                    if (field >= getMaxColumns())
                        setMaxColumns(field + 1);
                }
            } else {
                line[field] = value.toString();
            }
        }
    }

    /**
     * Sets the value of the given field as value.toString()
     * <br/>
     * If the row did not previously had enough columns to accommodate the field,
     * the needed amount of fields will be added with MISSING_VALUE.
     * 
     * @param value the new value
     * @param field column of the 
     */
    public void setValue(Object value, int field) {
        if (m_current >= m_data.size()) {
            throw new ArrayIndexOutOfBoundsException("Behind the last row, therefore I can't set a value");
        }
        if (m_current < 0) {
            throw new ArrayIndexOutOfBoundsException("Before the first row, therefore I can't set a value");
        }
        setValue(value, field, m_current);
    }
    
    
    /**
     * returns the specified value. 
     * If the field does not contain a value or the defined field is outside the CSV-file {@link #MISSING_FIELD} is returned
     * @param field column 
     * @param row row
     * @return value at the specified place.
     */
    public String getValue(Integer field, Integer row) {
        if (row == null || row >= m_data.size())
            return MISSING_FIELD;
        String[] line = m_data.get(row); 
        if (field == null || field >= line.length)
            return MISSING_FIELD;
        return line[field];
    }

    /**
     * returns the specified value interpreted as double. 
     * If the field does not contain a value, the defined field is outside the CSV-file or the value can not be interpreted as double Double.NaN is returned
     * @param field column 
     * @param row row
     * @return value at the specified place.
     */
    public double getDouble(Integer field, Integer row) {
        String v = getValue(field, row);
        if (v == MISSING_FIELD)
            return Double.NaN;
        try {
            return Double.parseDouble(v);
        } catch (NumberFormatException nfe) {
            return Double.NaN;
        }
    }
    
    /**
     * returns the specified value interpreted as boolean. 
     * <p>Any value that can not be matched by {@code ^(F|0|FALSE|N|NO)?$ } is taken as true.</p>
     * If the field does not contain a value or the defined field is outside the CSV-file false is returned
     * @param field column 
     * @param row row
     * @return value at the specified place.
     */
    public Boolean getBool(Integer field, Integer row) {
        String v = getValue(field, row);
        if (v == MISSING_FIELD)
            return false;
        
        //ISFALSE.split(v);
        return ! ISFALSE.matcher(v).matches();
    }

    /**
     * returns the specified value interpreted as boolean. 
     * <p>If the field does not contain a value or the defined field is outside the CSV-file the default value is returned.
     * </p> Otherwise the result depends on defaultValue.
     * <ul> 
     * <li>if defaultValue is true then anything that can not be matched by {@code ^(F|0|FALSE|N|NO|\-)?$ } is returned as true</li>
     * <li>if defaultValue is false only things that match {@code ^(T|1|-1|TRUE|Y|YES|\+)$ } is returned as true</li>
     * </ul>
     * @param field column 
     * @param row row
     * @param defaultValue the default value
     * @return value at the specified place.
     */
    public Boolean getBool(Integer field, Integer row, boolean defaultValue) {
        String v = getValue(field, row);
        if (v == MISSING_FIELD)
            return defaultValue;
        
        if (defaultValue)
            return !ISFALSE.matcher(v).matches();
        else
            return ISTRUE.matcher(v).matches();
    }

    
    /**
     * returns the specified value interpreted as double. 
     * If the field does not contain a value, the defined field is outside the CSV-file or the value can not be interpreted as double Double.NaN is returned
     * @param fieldName column 
     * @param row row
     * @return value at the specified place.
     */
    public double getDouble(String fieldName, Integer row) {
        Integer field = getColumn(fieldName);
        if (field == null)
            return Double.NaN;
        return getDouble(field, row);
    }
    
    /**
     * returns the specified value interpreted as Boolean. 
     * If the field does not contain a value, the defined field is outside the CSV-file or the value can not be interpreted as Boolean null is returned
     * @param fieldName column 
     * @param row row
     * @return value at the specified place.
     */
    public Boolean getBool(String fieldName, Integer row) {
        Integer field = getColumn(fieldName);
        if (field == null)
            return null;
        return getBool(field, row);
    }

    /**
     * returns the specified value interpreted as boolean. 
     * <p>If the field does not contain a value or the defined field is outside the CSV-file the default value is returned.
     * </p> Otherwise the result depends on defaultValue.
     * <ul> 
     * <li>if defaultValue is true then anything that can not be matched by {@code ^(F|0|FALSE|N|NO|\-)?$ } is returned as true</li>
     * <li>if defaultValue is false only things that match {@code ^(T|1|-1|TRUE|Y|YES|\+)$ } is returned as true</li>
     * </ul>
     * @param fieldName column 
     * @param row row
     * @param defaultValue the default value for miossing values
     * @return value at the specified place.
     */
    public Boolean getBool(String fieldName, Integer row, boolean defaultValue) {
        Integer field = getColumn(fieldName);
        return getBool(field, row,defaultValue);
    }

    /**
     * get the current row as a concatenated String
     * @return 
     */
    @Override
    public String getCurrentLine() {
        
        return getLine(m_current);
    }
    
    /**
     * get the defined row as a concatenated String
     * @param row the row to be return as string
     * @return 
     */
    public String getLine(int row) {
        String[] line = m_data.get(row);
        return valuesToString(line);        
    }
    
    /**
     * turns the given array into a CSV-line that fits to the currently defined delimiter and quote characters
     * @param values
     * @return 
     */
    public String valuesToString(String[] values) {
        StringBuilder sb = new StringBuilder(values.length*10);
        for (int i = 0; i< values.length; i++) {
            sb.append(quoteValue(values[i]));
            sb.append(getDelimiter());
        }
        return sb.substring(0, sb.length() -1);        
        
    }
    /**
     * turns the given array into a CSV-line that fits to the currently defined delimiter and quote characters
     * @param values
     * @return 
     */
    public String valuesToString(List<String> values) {
        String[] v =new String[values.size()];
        values.toArray(v);
        return valuesToString(v);
    }    
    /**
     * adds a new unnamed column to the csv-file
     * @return the index of the column
     */
    @Override
    public synchronized  int addColumn() {
        int ret = super.addColumn();
        notifyColumnsChanged();
        return ret;
    }

    /**
     * adds a new named column to the csv-file
     * @param name the name of the column
     * @return the index of the new column
     */
    @Override
    public synchronized int addColumn(String name) {
        int ret = super.addColumn(name);
        notifyColumnsChanged();
        return ret;
    }
    

    /**
     * Tests if a column of the given name exists and if it does it returns the index of the existing column.
     * If no column of the given name exists a new column is added and the index of the new column is returned.
     * @param name the name of the column
     * @return the index of the existing or newly created column
     */
    public synchronized int addGetColumn(String name) {
        Integer c = getColumn(name);
        if (c == null) {
            int ret = super.addColumn(name);
            notifyColumnsChanged();
            return ret;
        }
        return c;
    }


    /**
     * How many rows are in the file.
     * @return 
     */
    public int getRowCount() {
        return m_data.size();
    }
   
    /**
     * register a new listener, that is called when the file is finished reading in.
     * @param listener 
     */
    public void addListenerComplete(LoadListener listener) {
        if (!m_listenerCompleted.contains(listener))
            m_listenerCompleted.add(listener);
        else
            System.out.println("some error here");
    }
    
    /**
     * register a new listener, that is called regularly during loading of the file.
     * @param listener 
     */
    public void addListenerProgress(LoadListener listener) {
        m_listenerProgress.add(listener);
    }

    /**
     * removes the given listener
     * @param listener 
     */
    public void removeListenerComplete(LoadListener listener) {
        m_listenerCompleted.remove(listener);
    }

    /**
     * removes the given listener
     * @param listener 
     */
    public void removeListenerProgress(LoadListener listener) {
        m_listenerProgress.remove(listener);
    }

    /**
     * register a new listener, that is called when sorting is done
     * @param listener 
     */
    public void addListenerSort(LoadListener listener) {
        m_listenerSort.add(listener);
    }

    /**
     * register a new listener, that is called when columns where changed (added, removed, renamed)
     * @param listener 
     */
    public void addListenerColumnsChanged(LoadListener listener) {
        m_listenerColumnsChanged.add(listener);
    }

    /**
     * removes the given listener
     * @param listener 
     */
    public void removeListenerSort(LoadListener listener) {
        m_listenerSort.remove(listener);
    }

    /**
     * removes the given listener
     * @param listener 
     */
    public void removeListenerColumnsChanged(LoadListener listener) {
        m_listenerColumnsChanged.remove(listener);
    }

    /**
     * sort the file numeric by the given column
     * @param column 
     */
    public void sortNumeric(final int column) {
        if (column > getMaxColumns())
            return;
        Comparator<String[]> comp = new Comparator<String[]>(){
            
            public int compare(String[] s1, String[] s2) {
                if (s1.length <= column) {
                    if (s2.length <= column)
                        return 0;
                    else
                        return 1;
                }
                if (s2.length <= column)
                    return -1;
                
                Matcher m1 = NUMERIC.matcher(s1[column]);
                Matcher m2 = NUMERIC.matcher(s2[column]);
                if (!m1.matches()) {
                    if (!m2.matches())
                        return 0;
                    else
                        return 1;
                }
                if (!m2.matches())
                    return -1;
                double d1 = Double.parseDouble(m1.group(0));
                double d2 = Double.parseDouble(m2.group(0));
                    
                return Double.compare(d1, d2);
            }
        };
        synchronized(m_data) {
            java.util.Collections.sort(m_data, comp);
        }
        notifySort();
    }

    /**
     * sort the file numeric by the given column with the highest numbers first and the lowest last
     * @param column 
     */
    public void sortNumericReverse(final int column) {
        if (column > getMaxColumns())
            return;
        Comparator<String[]> comp = new Comparator<String[]>(){

            public int compare(String[] s1, String[] s2) {
                if (s1.length <= column) {
                    if (s2.length <= column)
                        return 0;
                    else
                        return -1;
                }
                if (s2.length <= column)
                    return 1;

                Matcher m1 = NUMERIC.matcher(s1[column]);
                Matcher m2 = NUMERIC.matcher(s2[column]);
                if (!m1.matches()) {
                    if (!m2.matches())
                        return 0;
                    else
                        return -1;
                }
                if (!m2.matches())
                    return 1;
                double d1 = Double.parseDouble(m1.group(0));
                double d2 = Double.parseDouble(m2.group(0));
                return Double.compare(d2, d1);
            }
        };
        synchronized(m_data) {
            java.util.Collections.sort(m_data, comp);
        }
        notifySort();
    }

    /**
     * sort the file according to the given sort-criteria
     * @param criteria 
     */
    public void sort(final CSVSort[] criteria) {
        Comparator<String[]> comp = new Comparator<String[]>(){
            
            public int compare(String[] s1, String[] s2) {
                for (int i =0; i< criteria.length;i++) {
                    int r = criteria[i].compare(s1, s2);
                    if (r!=0) 
                        return r;
                }
                return 0;
            }
        };
        synchronized(m_data) {
            java.util.Collections.sort(m_data, comp);
        }
        notifySort();
    }
    

    /**
     * sort the file revers to the given sort-criteria
     * @param criteria 
     */
    public void sortReverse(final CSVSort[] criteria) {
        Comparator<String[]> comp = new Comparator<String[]>(){
            
            public int compare(String[] s1, String[] s2) {
                for (int i =0; i< criteria.length;i++) {
                    int r = criteria[i].compare(s1, s2);
                    if (r!=0) 
                        return -r;
                }
                return 0;
            }
        };
        java.util.Collections.sort(m_data, comp);
        notifySort();
    }
    
    
    /**
     * sort the file alphanumerically by the given column
     * @param column column by which to sort the file 
     */
    public void sortAlpha(final int column) {
        if (column > getMaxColumns())
            return;
        Comparator<String[]> comp = new Comparator<String[]>(){
            
            public int compare(String[] s1, String[] s2) {
                if (s1.length <= column) {
                    if (s2.length <= column)
                        return 0;
                    else
                        return 1;
                }
                if (s2.length <= column)
                    return -1;
                
                return s1[column].compareTo(s2[column]) ;
            }
        };
        java.util.Collections.sort(m_data, comp);
        notifySort();
    }

    /**
     * sort the file reversely alphanumeric by the given column
     * @param column column by which to sort the file 
     */
    public void sortAlphaReverse(final int column) {
        if (column > getMaxColumns())
            return;
        Comparator<String[]> comp = new Comparator<String[]>(){
            
            public int compare(String[] s1, String[] s2) {
                if (s1.length <= column) {
                    if (s2.length <= column)
                        return 0;
                    else
                        return -1;
                }
                if (s2.length <= column)
                    return 1;
                
                return s2[column].compareTo(s1[column]) ;
            }
        };
        java.util.Collections.sort(m_data, comp);
        notifySort();
    }
    
    /**
     * make sure the given row has a value stored for each column. columns without values will get {@link #MISSING_FIELD} stored.
     * @param row 
     */
    public synchronized void ensureCapacity(int row) {
        String[] previous = m_data.get(row);
        int targetColumns = getMaxColumns();
        if (previous.length < targetColumns) {
            String[] ext = java.util.Arrays.copyOf(previous, targetColumns);
            java.util.Arrays.fill(ext, previous.length, targetColumns, MISSING_FIELD);
        }
    }

    /**
     * notify listeners that new data where read.
     * @param rows 
     */
    protected synchronized void  notifyProgress(int rows) {
        
        for (LoadListener l : (ArrayList<LoadListener>) m_listenerProgress.clone()) {
            l.listen(rows,-1);
        }
    }
    
    /**
     * notify listeners that all data where read.
     */
    protected synchronized void notifyComplete() {
        for (LoadListener l : (ArrayList<LoadListener>)m_listenerCompleted.clone()) {
            l.listen(m_data.size(),getMaxColumns());
        }
    }

    /**
     * notify listeners that one or more columns have changed
     */
    protected synchronized void notifyColumnsChanged() {
        for (LoadListener l : (ArrayList<LoadListener>)m_listenerColumnsChanged.clone()) {
            l.listen(m_data.size(),getMaxColumns());
        }
    }
    
    
    /**
     * notify listeners that sorting is finished
     */
    protected void notifySort() {
        synchronized(m_data) {
            for (LoadListener l : (ArrayList<LoadListener>)m_listenerSort.clone()) {
                l.listen(m_data.size(), -1);
            }
        }
    }

    /** 
     * 
     * @return if data are currently read.
     */
    public boolean isLoading() {
        return m_loading;
    }

    /**
     * write the data out to the given PrintWriter with the currently defined delimeter and quote characters.
     * @param out 
     */
    public void writeFile(PrintWriter out) {
        if (hasHeader()) {
            out.println(valuesToString(getHeader()));
        }
        for (String[] line : m_data) {
            out.println(valuesToString(line));
        }
    }

    /**
     * empty the data and reset the header
     */
    public void clear() {
        synchronized(m_data) {
            m_data.clear();
            setHeader("");
            notifyColumnsChanged();
        }
    }
    
    
    /**
     * if a header is defined, that header will be the first data-row and the header will be cleared
     * if no header is defined the first row will be turned into the header
     */
    public void switchHeader() {
        if (getHeader() == null) {
            if (!m_data.isEmpty()) {
                setHeader(m_data.get(0));
                deleteRow(0);
            }
        } else {
            insertLine(0, getHeader());
            cleanHeader();
        }
    }
    
    /**
     * add listeners to be notified if the header changes
     * @param l listener to add
     */
    public void addHeaderChangedListener(java.awt.event.ActionListener l) {
        m_listenerHeaderChanged.add(l);
    }
    
    /**
     * removes the given listener
     * @param l
     */
    public void removeHeaderChangedListener(java.awt.event.ActionListener l) {
        m_listenerHeaderChanged.remove(l);
    }
    
    /**
     * notify the listeners
     */
    protected void fireHeaderChanged() {
        ActionEvent evt = new ActionEvent(this,ActionEvent.ACTION_PERFORMED,"HEADER CHANGED" );
        for (java.awt.event.ActionListener e : m_listenerHeaderChanged) {
            e.actionPerformed(evt);
        }
                    
    }

    /**
     * set the column-headers according to the array
     * @param header 
     */
    public void setHeader(String[] header) {
        super.setHeader(header);
        fireHeaderChanged();
    }
    
    /** clear headers*/
    public void cleanHeader() {
        super.cleanHeader();
        fireHeaderChanged();
    }
    
}
