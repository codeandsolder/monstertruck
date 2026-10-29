//! Converts OBJ and STL to each other.
//!
//! usage:
//!
//! ```bash
//! cargo run --example obj_stl <input-file>
//! ```

use monstertruck_mesh::{
    obj,
    stl::{self, StlType},
};
use std::{error::Error, fs::File, io, path::PathBuf};

fn main() -> Result<(), Box<dyn Error>> {
    let path = std::env::args_os()
        .nth(1)
        .map(PathBuf::from)
        .ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidInput, "usage: obj_stl <input-file>")
        })?;
    let extension = path.extension().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("cannot infer file type from {}", path.display()),
        )
    })?;

    if extension == "obj" {
        let polymesh = obj::read(File::open(&path)?)?;
        let output = path.with_extension("stl");
        let mut file = File::create(output)?;
        stl::write(&polymesh, &mut file, StlType::Binary)?;
    } else if extension == "stl" {
        let polymesh = stl::read(File::open(&path)?, StlType::Automatic)?;
        let output = path.with_extension("obj");
        let mut file = File::create(output)?;
        obj::write(&polymesh, &mut file)?;
    } else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("unsupported file type for {}", path.display()),
        )
        .into());
    }

    Ok(())
}
