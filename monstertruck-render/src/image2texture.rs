use crate::{
    DeviceHandler, Extent3d, Origin3d, RenderError, TexelCopyBufferLayout, TexelCopyTextureInfo,
    Texture, TextureAspect, TextureDescriptor, TextureDimension, TextureFormat, TextureUsages,
};
use image::{DynamicImage, RgbaImage};

/// Utility for creating `Texture` from `DynamicImage`
#[inline]
/// # Errors
/// Returns [`RenderError::InvalidTextureDimensions`] when the image has a zero-sized axis or exceeds the GPU's 2D texture limit, and [`RenderError::TextureRowTooLarge`] when the row stride cannot fit wgpu's `u32` layout field.
pub fn image2texture(
    device_handler: &DeviceHandler,
    image: &DynamicImage,
) -> Result<Texture, RenderError> {
    let buffer = image.to_rgba8();
    imagebuffer2texture(device_handler, &buffer)
}

const fn validate_texture_dimensions(
    width: u32,
    height: u32,
    max_dimension: u32,
) -> Result<(), RenderError> {
    if width == 0 || height == 0 || width > max_dimension || height > max_dimension {
        Err(RenderError::InvalidTextureDimensions {
            width,
            height,
            max_dimension,
        })
    } else {
        Ok(())
    }
}

fn texture_bytes_per_row(width: u32, pixel_size: usize) -> Result<u32, RenderError> {
    let pixel_size_u32 = u32::try_from(pixel_size)
        .map_err(|_| RenderError::TextureRowTooLarge { width, pixel_size })?;
    width
        .checked_mul(pixel_size_u32)
        .ok_or(RenderError::TextureRowTooLarge { width, pixel_size })
}

fn imagebuffer2texture(
    device_handler: &DeviceHandler,
    image_buffer: &RgbaImage,
) -> Result<Texture, RenderError> {
    let (device, queue) = (device_handler.device(), device_handler.queue());
    let size = Extent3d {
        width: image_buffer.width(),
        height: image_buffer.height(),
        depth_or_array_layers: 1,
    };
    validate_texture_dimensions(
        size.width,
        size.height,
        device.limits().max_texture_dimension_2d,
    )?;
    let texture = device.create_texture(&TextureDescriptor {
        label: None,
        size,
        mip_level_count: 1,
        sample_count: 1,
        dimension: TextureDimension::D2,
        format: TextureFormat::Rgba8Unorm,
        view_formats: &[],
        usage: TextureUsages::TEXTURE_BINDING | TextureUsages::COPY_DST,
    });
    let bytes_per_row = texture_bytes_per_row(size.width, 4)?;
    queue.write_texture(
        TexelCopyTextureInfo {
            texture: &texture,
            mip_level: 0,
            origin: Origin3d::ZERO,
            aspect: TextureAspect::All,
        },
        image_buffer.as_raw(),
        TexelCopyBufferLayout {
            offset: 0,
            bytes_per_row: Some(bytes_per_row),
            rows_per_image: Some(size.height),
        },
        size,
    );
    Ok(texture)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn texture_dimensions_reject_zero_axes() {
        assert!(validate_texture_dimensions(0, 1, 8192).is_err());
        assert!(validate_texture_dimensions(1, 0, 8192).is_err());
    }

    #[test]
    fn texture_dimensions_reject_device_limit_overflow() {
        assert!(validate_texture_dimensions(8193, 1, 8192).is_err());
        assert!(validate_texture_dimensions(1, 8193, 8192).is_err());
    }

    #[test]
    fn texture_dimensions_accept_device_limit_boundary() {
        assert!(validate_texture_dimensions(8192, 8192, 8192).is_ok());
    }

    #[test]
    fn row_stride_rejects_u32_multiplication_overflow() {
        assert!(texture_bytes_per_row(u32::MAX, 2).is_err());
    }

    #[cfg(target_pointer_width = "64")]
    #[test]
    fn row_stride_rejects_pixel_sizes_wider_than_u32() {
        let pixel_size = 1usize << 32;
        assert!(texture_bytes_per_row(1, pixel_size).is_err());
    }

    #[test]
    fn row_stride_accepts_representable_layouts() {
        assert!(matches!(texture_bytes_per_row(1920, 4), Ok(7680)));
    }
}
