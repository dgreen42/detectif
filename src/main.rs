// In this I will attempt to write a program tht can identify a pattern that shows up more than
// once in a specific sequence

use core::str;
use std::{
    collections::HashMap,
    env::args,
    fs::{read, File},
    io::{BufRead, Write},
};

fn main() {
    let input = args().nth(1).expect("Enter a file to check");
    let file = read(input).unwrap();
    for line in file.lines() {
        search(line.unwrap(), 8);
    }
}

fn search(sequence: String, window_cap: u16) -> HashMap<String, u32> {
    let map: HashMap<String, u32> = HashMap::new();
    // I don't think we need ot specify what strand we are on since this will be lookig at raw
    // nano pore sequence. Strand would only matter for mapping

    let bstring = sequence.as_bytes();
    let window_size: u16 = 3;

    while window_size < window_cap {
        for begin in 0..(bstring.len() - (window_size as usize)) {
            let check = map_check(&bstring[begin..(window_size as usize)], &map);
        }
    }

    map
}

fn map_check(pattern: &[u8], hashmap: &HashMap<String, u32>) -> u8 {
    let keys = hashmap.keys();
    let mut matched: u8 = 0;
    for key in keys {
        if key == str::from_utf8(pattern).unwrap() {
            matched += 1;
        }
    }
    matched
}
