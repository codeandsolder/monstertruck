//! Per-surface resident GPU seed search. Optional and strictly cross-checked
//! against the FP64 CPU index before its result can reach the CAD solver.

use super::{Bounded, MetricSpace};
use std::{ffi::c_void, ptr::NonNull};

unsafe extern "C" {
    fn mt_gpu_grid_create(
        positions: *const f64,
        count: u32,
        dimensions: u32,
        output: *mut *mut c_void,
    ) -> i32;
    fn mt_gpu_grid_nearest(handle: *mut c_void, target: *const f64, result: *mut u32) -> i32;
    fn mt_gpu_grid_release(handle: *mut c_void);
}

pub(super) struct CudaSeedGrid {
    handle: NonNull<c_void>,
}

impl CudaSeedGrid {
    pub(super) fn new<P>(samples: &[(f64, f64, P)]) -> Option<Self>
    where P: Bounded<Scalar = f64> + MetricSpace<Metric = f64> + Copy {
        let dimensions = samples.first()?.2.dimension();
        if dimensions == 0 || dimensions > 4 || samples.len() > u32::MAX as usize {
            return None;
        }
        let mut positions = Vec::with_capacity(samples.len() * 4);
        for (_, _, point) in samples {
            for axis in 0..4 {
                positions.push(if axis < dimensions { point[axis] } else { 0.0 });
            }
        }
        if positions.iter().any(|value| !value.is_finite()) {
            return None;
        }
        let mut raw = std::ptr::null_mut();
        // SAFETY: positions is contiguous and remains borrowed for the entire
        // synchronous CUDA grid-upload call; output points to a valid pointer.
        let status = unsafe {
            mt_gpu_grid_create(
                positions.as_ptr(),
                u32::try_from(samples.len()).ok()?,
                u32::try_from(dimensions).ok()?,
                &mut raw,
            )
        };
        if status == 0 {
            Some(Self {
                handle: NonNull::new(raw)?,
            })
        } else {
            None
        }
    }

    pub(super) fn nearest<P>(&self, target: P) -> Option<usize>
    where P: Bounded<Scalar = f64> + MetricSpace<Metric = f64> + Copy {
        let dimensions = target.dimension();
        let mut point = [0.0_f64; 4];
        if dimensions == 0 || dimensions > 4 {
            return None;
        }
        for (axis, entry) in point.iter_mut().enumerate().take(dimensions) {
            *entry = target[axis];
        }
        if point.iter().any(|value| !value.is_finite()) {
            return None;
        }
        let mut result = u32::MAX;
        // SAFETY: handle owns a live CUDA grid, this thread owns the source
        // and result buffers, and the native call synchronizes before return.
        let status =
            unsafe { mt_gpu_grid_nearest(self.handle.as_ptr(), point.as_ptr(), &mut result) };
        (status == 0 && result != u32::MAX).then_some(result as usize)
    }
}

impl Drop for CudaSeedGrid {
    fn drop(&mut self) {
        // SAFETY: handle was returned by mt_gpu_grid_create and is freed once.
        unsafe { mt_gpu_grid_release(self.handle.as_ptr()) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::Point3;

    #[test]
    fn gpu_grid_matches_ordered_fp64_distance_scan() {
        let samples = (0..=50)
            .flat_map(|u| {
                (0..=50).map(move |v| {
                    (
                        0.0,
                        0.0,
                        Point3::new(
                            67.0 + f64::from(u) / 50.0,
                            5.0 + f64::from(v) / 50.0,
                            1.7 + (f64::from(v) / 50.0).sin(),
                        ),
                    )
                })
            })
            .collect::<Vec<_>>();
        let gpu = CudaSeedGrid::new(&samples).expect("GPU should accept f64 sample grid");
        for i in (0..samples.len()).step_by(7) {
            let probe = samples[i].2;
            let cpu = samples
                .iter()
                .enumerate()
                .fold((f64::INFINITY, None), |best, (idx, p)| {
                    let d = p.2.distance2(probe);
                    if d < best.0 { (d, Some(idx)) } else { best }
                })
                .1;
            assert_eq!(gpu.nearest(probe), cpu);
        }
    }
}
