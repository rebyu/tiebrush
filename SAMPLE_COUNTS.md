## Running TieBrush with exact sample counts

This fork of TieBrush includes updates to be able to track exact sample counts when using the flag `-C` or `--sample-counts`.

### New structs/classes created in order to track exact sample counts

- SCounts (struct)

    _parameters_
    ```c++
    int tid                 // chromosome key #
    std::string chr_name    // chromosome name (retrieved from SAM header given tid)
    int start               // start location of sample
    int end                 // end location of sample
    int scount              // number of samples from start-end
    bool from_tb            // is this sample retrieved from 
    std::set<int> samp_ids  // set of samples including a read from start-end
    ```

- SCData (class)

    _parameters_
    ```c++
    std::vector<SCounts*> sclst     // vector of SCounts objects, buffer as samples are read in
    ```
    _functions_
    ```c++
    void clear()                    // function to clear sclst
    void collapse()                 // collapse sclst to contain accumulated sample counts
    void mergeCounts(int tid)       // merge existing sclst counts with counts from tiebrush output files given a tid
    void add(GList<SPData>& spdlst) // add data from spdlst to sclst as SCounts* object
    void flush(int outtid)          // output to file and flush sclst given a tid
    ```


### General Logic

- Add `SCounts` object(s) to vector given updated spdlst every time the input reaches a new position (assuming non-decreasing input)
- Flush `sclst` (vector of `SCounts`) every time a new chromosome is reached and at end of parsing
- All TieBrush output files have sample counts loaded into a map of `<"chromosome_name":std::vector<SCounts*>>` that are compiled into a vector (`tbcounts`)


Overall, to track exact sample counts, modifications were made to `tiebrush.cpp`, `tmerge.cpp`, `tmerge.h`.

- Minor modifications in `tmerge` were to track file names while loading samples in for tracking sample counts of TieBrush output files.

- `tiewrap.py` was also updated to take in `-C` and `--sample-counts` as a parameter in order for `run_tiebrush.py` to batch run TieBrush with exact sample count tracking.