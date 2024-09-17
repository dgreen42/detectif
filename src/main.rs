// In this I will attempt to write a program tht can identify a pattern that shows up more than
// once in a specific sequence

use core::str;
use csv::Writer;
use std::{
    collections::HashMap,
    env::args,
    fs::{read, File},
    io::{BufRead, Write},
    thread::current,
    usize,
};

fn main() {
    let input = args().nth(1).expect("Enter a file to check");
    let path = args().nth(2).expect("Enter output file name");
    let file = read(input).unwrap();
    let mut aggregate_lines = String::new();
    for line in file.lines() {
        aggregate_lines.push_str(&line.unwrap());
    }
    let map = search(aggregate_lines, 8);
    write_csv(&map, path);
}

fn write_csv(hashmap: &HashMap<String, u32>, path: String) {
    let mut writer = Writer::from_path(path).expect("Could not open file");
    for key in hashmap.keys() {
        let current_pair = hashmap.get_key_value(key).unwrap();
        let mut current_record: Vec<String> = Vec::new();
        current_record.push(current_pair.0.to_string());
        current_record.push(current_pair.1.to_string());
        writer.write_record(current_record).unwrap();
    }
}

fn search(sequence: String, window_cap: u16) -> HashMap<String, u32> {
    let mut map: HashMap<String, u32> = HashMap::new();
    // I don't think we need ot specify what strand we are on since this will be lookig at raw
    // nano pore sequence. Strand would only matter for mapping

    let bstring = sequence.as_bytes();
    let mut window_size: u16 = 3;

    while window_size < window_cap {
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
