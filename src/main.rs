use ::image as im;
use piston_window::*;
use cgmath::{dot, Vector2, Vector3, Matrix2, Rad, Deg, InnerSpace};
use fps_counter::FPSCounter;
use rand::Rng;

const WINDOW_WIDTH: u32 = 800;
const WINDOW_HEIGHT: u32 = 400;

const WINDOW_FWIDTH: f64 = WINDOW_WIDTH as f64;
const WINDOW_FHEIGTH: f64 = WINDOW_HEIGHT as f64;

const PI: f64 = 3.1415926535897932385;

fn random_0_1() -> f64 {
    return rand::thread_rng().gen_range(0.0..1.0)
}

fn random_in_interval(i: Interval) -> f64 {
    return rand::thread_rng().gen_range(i.min..i.max)
}

struct Interval {
    min: f64,
    max: f64
}

impl Interval {
    fn new() -> Interval {
        Interval {min: f64::INFINITY, max: f64::NEG_INFINITY}
    }
    fn from(min: f64, max: f64) -> Interval{
        Interval {min: min, max: max}
    }

    fn contains(&self, n: f64) -> bool {
        return self.min <= n && n <= self.max
    }

    fn surrounds(&self, n: f64) -> bool {
        return self.min < n && n < self.max
    }

    fn clamp(&self, n: f64) -> f64 {
        if n < self.min {return self.min}
        if n > self.max {return self.max}
        n
    }
}

const EMPTY: Interval = Interval {min: f64::INFINITY, max: f64::NEG_INFINITY};

const UNIVERSE: Interval = Interval {min: f64::NEG_INFINITY, max: f64::INFINITY};


#[derive(Clone, Copy)]
struct HitRecord {
    position: Vector3<f64>,
    normal: Vector3<f64>,
    t: f64,
    front_face: bool
}

impl HitRecord {
    fn new() -> HitRecord {
        HitRecord {position: Vector3::from([0.,0.,0.]), normal: Vector3::from([0.,0.,0.]), t: 0., front_face: true}
    }

    fn set_face_normal(&mut self, r: Ray, outward_normal: Vector3<f64>) {
        self.front_face = dot(r.direction, outward_normal) < 0.;
        self.normal = if self.front_face {outward_normal} else {-outward_normal};
    }
}
#[derive(Clone)]
struct Sphere {
    position: Vector3<f64>,
    radius: f64
}
#[derive(Clone)]
struct World {
    spheres: Vec<Sphere>
}

impl World {
    fn new() -> World {
        World {spheres: vec![]}
    }

    fn hit_spheres(&self, r: Ray, ray_t: Interval, rec: &mut HitRecord) -> bool {
        let mut temp_rec = HitRecord::new();
        let mut hit_anything: bool = false;
        let mut closest_so_far = ray_t.max;

        for sphere in &self.spheres {
            if r.hit_sphere(sphere.position, sphere.radius, Interval::from(ray_t.min, closest_so_far), &mut temp_rec) {
                hit_anything = true;
                closest_so_far = temp_rec.t;
                *rec = temp_rec;
            }
        }
        hit_anything
    }
}

#[derive(Clone, Copy)]
struct Ray {
    origin: Vector3<f64>,
    direction: Vector3<f64>
}

impl Ray {
   fn new(origin: Vector3<f64>, direction: Vector3<f64>) -> Ray{
        Ray {origin, direction}
   }
   fn at(&self, t: f64) -> Vector3<f64> {
        self.origin + t * self.direction
   }
   fn hit_sphere(&self, sphere_center: Vector3<f64>, sphere_radius: f64, ray_t: Interval, rec: &mut HitRecord) -> bool {
        let oc = self.origin - sphere_center;
        let a: f64 = self.direction.magnitude2();
        let half_b = dot(oc, self.direction);
        let c = oc.magnitude2() - sphere_radius * sphere_radius;
        let discriminant = half_b*half_b -a*c;
        if discriminant < 0. {
            return false
        } 
        let sqrtd = discriminant.sqrt();
        let mut root = (-half_b - sqrtd) / a;
        if !ray_t.surrounds(root) {
            root = (-half_b + sqrtd) / a;
            if !ray_t.surrounds(root) { 
                return false
            }
        }
        rec.t = root;
        rec.position = self.at(rec.t);
        rec.normal = (rec.position - sphere_center) / sphere_radius;
        let outward_normal = (rec.position - sphere_center) / sphere_radius;
        rec.set_face_normal(*self, outward_normal);

        true

   }

   fn ray_color(&self, world: World) -> [u8; 4] {
        let mut rec = HitRecord::new();
        if world.hit_spheres(*self, Interval::from(0., f64::INFINITY) , &mut rec) {
            let r = (0.5 * (rec.normal.x+1.))*255.99;
            let g = (0.5 * (rec.normal.y+1.))*255.99;
            let b = (0.5 * (rec.normal.z+1.))*255.99;
            return [r as u8, g as u8, b as u8, 255]
        }
        let a = (self.direction.normalize().y + 1.) / 2.;
        let r = ((1.-a) * 1. + a * 0.5)*255.99;
        let g = ((1.-a) * 1. + a * 0.7)*255.99;
        let b = ((1.-a) * 1. + a * 1.)*255.99;
        [r as u8, g as u8, b as u8, 255]


   }
   //fn color(&self) -> [u8; 4] {
   //     let t = self.hit_sphere(Vector3::from([0., 0., -1.]), 0.5);
   //     if t > 0. {
   //         let normal = (self.at(t) - Vector3::from([0., 0., -1.])).normalize();
   //         let r = (0.5 * (normal.x+1.))*255.99;
   //         let g = (0.5 * (normal.y+1.))*255.99;
   //         let b = (0.5 * (normal.z+1.))*255.99;
   //         return [r as u8, g as u8, b as u8, 255]
   //         
   //     }
   //     let a = (self.direction.normalize().y + 1.) / 2.;
   //     let r = ((1.-a) * 1. + a * 0.5)*255.99;
   //     let g = ((1.-a) * 1. + a * 0.7)*255.99;
   //     let b = ((1.-a) * 1. + a * 1.)*255.99;
   //     [r as u8, g as u8, b as u8, 255]
   //}
}

fn render(canvas: &mut im::ImageBuffer<im::Rgba<u8>, Vec<u8>>, camera_position: Vector3<f64>, world: World) {

    let focal_length: f64 = 1.;
    let viewport_width = 4.;
    let viewport_height: f64 = 2.;

    let viewport_u: Vector3<f64> = Vector3::from([viewport_width, 0., 0.]);
    let viewport_v: Vector3<f64> = Vector3::from([0., -viewport_height, 0.]);

    let pixel_delta_u = viewport_u / WINDOW_FWIDTH;
    let pixel_delta_v = viewport_v / WINDOW_FHEIGTH;

    let viewport_top_left = camera_position - Vector3::from([0., 0., focal_length]) - viewport_u/2. - viewport_v/2.;
    let pixel00_loc = viewport_top_left + (pixel_delta_u + pixel_delta_v) / 2.;


    for x in 0..WINDOW_WIDTH {
        for y in 0..WINDOW_HEIGHT {
            let pixel_center = pixel00_loc + pixel_delta_u * x as f64 + pixel_delta_v * y as f64;
            let ray_direction = pixel_center - camera_position;

            let ray = Ray::new(camera_position, ray_direction);
            let r = x as f32 / (WINDOW_WIDTH as f32-1.)*255.99;
            let g = y as f32 / (WINDOW_HEIGHT as f32-1.)*255.99;
            let b = 0;
            canvas.put_pixel(x, y, im::Rgba(ray.ray_color(world.clone())));
        }
    }          
}


fn main() {
    let opengl = OpenGL::V3_2;
    let mut window: PistonWindow =
        WindowSettings::new("", (WINDOW_WIDTH, WINDOW_HEIGHT))
        .exit_on_esc(true)
        .graphics_api(opengl)
        .build()
        .unwrap();

    let mut fps_counter = FPSCounter::new();

    let mut canvas = im::ImageBuffer::new(WINDOW_WIDTH, WINDOW_HEIGHT);
    let mut texture_context = TextureContext {
        factory: window.factory.clone(),
        encoder: window.factory.create_command_buffer().into()
    };
    let mut texture: G2dTexture = Texture::from_image(
            &mut texture_context,
            &canvas,
            &TextureSettings::new()
        ).unwrap();


    let camera_position: Vector3<f64> = Vector3::from([0., 0., 0.]);
    let mut world = World::new();
    world.spheres.push(Sphere { position: Vector3::from([0., 0., -1.]), radius: 0.5 });
    world.spheres.push(Sphere { position: Vector3::from([0., -100.5, -1.]), radius: 100. });

    while let Some(e) = window.next() {
        let fps = fps_counter.tick() as f64;
        dbg!(fps);
        if e.render_args().is_some() {
            texture.update(&mut texture_context, &canvas).unwrap();
            window.draw_2d(&e, |c, g, device| {
                // Update texture before rendering.
                texture_context.encoder.flush(device);

                clear([0.0; 4], g);

                image(&texture, c.transform, g);
            });
        }

        canvas = im::ImageBuffer::new(WINDOW_WIDTH, WINDOW_HEIGHT);
        render(&mut canvas, camera_position, world.clone());
    }
}

