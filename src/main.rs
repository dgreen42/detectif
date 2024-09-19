// In this I will attempt to write a program tht can identify a pattern that shows up more than
// once in a specific sequence

use core::str;
use csv::Writer;
use indicatif::ProgressBar;
use std::{
    collections::HashMap,
    env::{args, set_var},
    fs::read,
    io::{stdin, BufRead},
    time::Duration,
    usize,
};

fn main() {
    set_var("RUST_BACKTRACE", "1");
    let in_opt = args().nth(1).expect("Enter the input option");
    if in_opt == "--help" {
        println!(
            "
detectif v0.0.1

detectif [optiont] [path] [window cap] [output.csv]

Example: detecif -f path/to/file 18 output.csv

Goes through a sequence and finds patterns. Made to find motifs in a nucleotide sequence. Starts at a base length of 3 nucleotides and increases the size of the window until it reaches a user defined cap. 

By itself detectif takes a fasta file and searches through the sequence outputtin a .csv with the resulting patterns and the number of times that they show up. 

In the genome pipeline it can take the sequences as a standard input. Using mokuba with the -sio flag, sequnces can be piped in. The only problem at the moment is that, if there are several input sequences then you should get multiple out. Still need to write this, but it is one the way.

"
        );
        std::process::exit(3);
    }
    let mut window_cap = String::new();
    let mut path = String::new();
    if in_opt == "-si" {
        let mut input: Vec<String> = Vec::new();
        let mut ids: Vec<String> = Vec::new();
        let std_input = stdin().lines();
        for line in std_input {
            let line = line.unwrap();
            if line.starts_with(">") {
                ids.push(line[1..].to_string());
            } else {
                input.push(line);
            }
        }
        window_cap.push_str(&args().nth(2).expect("Enter output Maximum motif length"));
        path.push_str(&args().nth(3).expect("Enter output file name"));
        let mut count = 1;
        let total = input.clone();
        for id in input {
            let map = search(
                &id,
                window_cap.parse::<u16>().unwrap(),
                &ids[count - 1],
                total.len() as i32,
                count as i32,
            );
            let numbered_path = format!("{}{}.csv", path, &ids[count - 1]);
            write_csv(&map, &numbered_path);
            count += 1;
        }
    }

    if in_opt == "-f" {
        let mut input = String::new();
        input.push_str(&args().nth(2).expect("Enter a file to check"));
        window_cap.push_str(&args().nth(3).expect("Enter output Maximum motif length"));
        path.push_str(&args().nth(4).expect("Enter output file name"));
        path.push_str(".csv");
        let file = read(input).unwrap();
        let mut aggregate_lines = String::new();
        let mut id = String::new();
        for line in file.lines() {
            let line = line.unwrap();
            if line.starts_with(">") {
                id.push_str(&line);
            }
            aggregate_lines.push_str(&line);
        }
        let map = search(
            &aggregate_lines,
            window_cap.parse::<u16>().unwrap(),
            &id,
            1,
            1,
        );
        write_csv(&map, &path);
    }
}

fn write_csv(hashmap: &HashMap<String, u32>, path: &str) {
    let mut writer = Writer::from_path(path).expect("Could not open file");
    writer
        .write_record(Vec::from(["pattern", "count"]))
        .unwrap();
    for key in hashmap.keys() {
        let current_pair = hashmap.get_key_value(key).unwrap();
        let mut current_record: Vec<String> = Vec::new();
        current_record.push(current_pair.0.to_string());
        current_record.push(current_pair.1.to_string());
        writer.write_record(current_record).unwrap();
    }
}

fn search(
    sequence: &str,
    window_cap: u16,
    id: &str,
    total: i32,
    count: i32,
) -> HashMap<String, u32> {
    let mut map: HashMap<String, u32> = HashMap::new();
    // I don't think we need ot specify what strand we are on since this will be lookig at raw
    // nano pore sequence. Strand would only matter for mapping

    let lowercase = sequence.to_lowercase();
    let bstring = lowercase.as_bytes();
    let mut window_size: u16 = 3;

    let spinner = ProgressBar::new_spinner()
        .with_message(format!("[{}/{}]Finding patters in {}", count, total, id));
    spinner.enable_steady_tick(Duration::from_millis(100));
    while window_size < window_cap {
        if sequence == "NA" {
            break;
        }
        for begin in 0..(bstring.len() - (window_size as usize)) {
            let current_pattern = str::from_utf8(&bstring[begin..begin + window_size as usize])
                .unwrap()
                .to_string();
            let check = map_check(&current_pattern, &map);
            if check == 0 {
                map.insert(current_pattern, 1);
            } else {
                let kv = map.get_key_value(&current_pattern).unwrap();
                let value = kv.1 + 1;
                map.insert(current_pattern, value);
            }
        }
        window_size += 1;
    }
    spinner.finish();
    map
}

fn map_check(pattern: &str, hashmap: &HashMap<String, u32>) -> u8 {
    let keys = hashmap.keys();
    let mut matched: u8 = 0;
    for key in keys {
        if key == pattern {
            matched += 1;
        }
    }
    matched
}
