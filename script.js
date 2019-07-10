let canvas = document.querySelector('#render');
let ctx = canvas.getContext('2d');

function resize() {
    canvas.style.width = window.innerWidth + 'px';
    canvas.style.height = window.innerHeight + 'px';
    canvas.width = canvas.clientWidth;
    canvas.height = canvas.clientHeight;
}
window.addEventListener('resize', resize);
resize();

class Vec2 {
    constructor(x, y) {
        this.x = x;
        this.y = y;
    }

    set(x, y) {
        this.x = x;
        this.y = y;
    }

    set(v) {
        this.x = v.x;
        this.y = v.y;
    }

    addVec(v) {
        this.x += v.x;
        this.y += v.y;
    }

    addMult(v, s) {
        this.addVec(this.multiplied(v, s));
    }

    lengthSq() {
        return this.x * this.x + this.y * this.y;
    }

    static multiplied(v, s) {
        return new Vec2(v.x * s, v.y * s);
    }

    static cross(a, b) {
        return a.x * b.y - a.y * b.x;
    }

    static newArray(length) {
        let arr = [];
        for (let i = 0; i < length; ++i) {
            arr.push(new Vec2(0, 0));
        }
        return arr;
    }
}

class Mat2 {
    constructor() {
        this.m00 = 0;
        this.m01 = 0;
        this.m10 = 0;
        this.m11 = 0;
    }

    set(radians) {
        let c = Math.cos(radians);
        let s = Math.sin(radians);

        this.m00 = c;
        this.m01 = -s;
        this.m10 = s;
        this.m11 = c;
    }

    set(mat) {
        this.m00 = mat.m00;
        this.m01 = mat.m01;
        this.m10 = mat.m10;
        this.m11 = mat.m11;
    }
}

const ImpulseMath = {
    EPSILON: 0.0001,
    EPSILON_SQ: () => {
        return this.EPSILON * this.EPSILON;
    },
    BIAS_RELATIVE: 0.95,
    BIAS_ABSOLUTE: 0.01,
    DT: 1 / 60,
    GRAVITY: new Vec2(0, 50),
    RESTING: () => {
        return multiplied(this.GRAVITY, this.DT).lengthSq() + this.EPSILON;
    },
    PENETRATION_ALLOWANCE: 0.05,
    PENETRATION_CORRECTION: 0.4,
    equal: (a, b) => {
        return Math.abs(a - b) <= this.EPSILON;
    },
    clamp: (min, max, val) => {
        return (val < min ? min : (val > max ? max : val));
    },
    random: (min, max) => {
        return (max - min) * Math.random() + min;
    },
    gt: (a, b) => {
        return a >= (b * this.BIAS_RELATIVE + a * this.BIAS_ABSOLUTE);
    }
};

class Shape {
    constructor(params) {
        this.type = params.type;
        this.radius = 0;
        this.u = new Mat2();
        // Obviously there's a better way to do this in JS (using subclasses and such),
        // but I'm trying to keep it as Scratch-like as possible to make it easier to port.
        if (this.type === 'c') {
            this.radius = params.radius;
        }
        else if (this.type === 'p') {
            this.vertices = Vec2.newArray(8);
            this.normals = Vec2.newArray(8);
            this.vertexCount = 0;
            if (params.hw) {
                vertexCount = 4;
                vertices[0].set(-params.hw, -params.hh);
                vertices[1].set(params.hw, -params.hh);
                vertices[2].set(params.hw, params.hh);
                vertices[3].set(-params.hw, params.hh);
                normals[0].set(0, -1);
                normals[1].set(1, 0);
                normals[2].set(0, 1);
                normals[3].set(-1, 0);
            }
            /*else {
                let verts = params.vertices;
                // Find the right most point on the hull
                let rightMost = 0;
                let highestX = verts[0].x;
                for (let i = 1; i < verts.length; ++i) {
                    let x = verts[i].x;
                    if (x > highestX) {
                        highestX = x;
                        rightMost = i;
                    }
                    // If same x, take lowest y
                    else if (x == highestX) {
                        if (verts[i].y < verts[rightMost].y) {
                            rightMost = i;
                        }
                    }
                }

                let hull = new Array(8);
                let outCount = 0;
                let indexHull = rightMost;

                for (;;) {
                    hull[outCount] = indexHull;

                    // Find next index that wraps around hull (find the next most counter-clockwise vertex)
                    let nextHullIndex = 0;
                    for (let i = 1; i < verts.length; ++i) {
                        // Skip if the same coord
                        if (nextHullIndex = indexHull) {
                            nextHullIndex = i;
                            continue;
                        }

                        // Cross every set of three unique verts
                    }
                }
            }*/
        }
    }

    clone() {
        if (this.type === 'c') {
            return new Shape({
                type: 'c',
                radius: this.radius
            });
        }
        else if (this.type === 'p') {
            let s = new Shape({
                type: 'p'
            });
            s.u.set(this.u);
            for (let i = 0; i < this.vertexCount; ++i) {
                s.vertices[i].set(this.vertices[i]);
                s.normals[i].set(this.normals[i]);
            }
            s.vertexCount = this.vertexCount;
            return s;
        }
    }

    initialize() {
        this.computeMass(1);
    }

    computeMass(density) {
        if (this.type === 'c') {
            this.body.mass = Math.PI * this.radius * this.radius * density;
            this.body.invMass = (this.body.mass !== 0) ? (1 / this.body.mass) : 0;
            this.body.inertia = this.body.mass * this.radius * this.radius;
            this.body.invInertia = (this.body.inertia !== 0) ? (1 / this.body.inertia) : 0;
        }
        else if (this.type === 'p') {
            //We'll get there...
        }
    }
    
    setOrient(radians) {
        if (this.type === 'p') {
            this.u.set(radians);
        }
    }
    
    /*getType() {
        return this.type;
    }*/
}

class Body {
    constructor(shape, x, y) {
        this.position = new Vec2(x, y);
        this.velocity = new Vec2(0, 0);
        this.force = new Vec2(0, 0);
        this.orient = (Math.random() - 0.5) * Math.PI;
        this.angularVelocity = 0;
        this.torque = 0;
        this.staticFriction = 0.5;
        this.dynamicFriction = 0.3;
        this.restitution = 0.2;
        this.mass = 0;
        this.invMass = 0;
        this.inertia = 0;
        this.invInertia = 0;
        this.shape = shape;
        shape.body = this;
        shape.initialize();
    }

    applyForce(force) {
        this.force.addVec(force);
    }

    applyImpulse(impulse, contactVector) {
        this.velocity.addMult(impulse, this.invMass);
        this.angularVelocity += this.invInertia * Vec2.cross(contactVector, impulse);
    }

    setStatic() {
        this.mass = 0;
        this.invMass = 0;
        this.inertia = 0;
        this.invInertia = 0;
    }

    setOrient(radians) {
        orient = radians;
        this.shape.setOrient(radians);
    }

    renderBody(ctx) {
        if (this.shape.type === 'c') { //Again, this isn't the best way, but I'm keeping it closer to Scratch
            ctx.beginPath();
            ctx.arc(this.position.x, this.position.y, this.shape.radius, 0, 2 * Math.PI);
            ctx.stroke();
        }
    }
}

class ImpulseScene {
    constructor(dt, iterations) {
        this.dt = dt;
        this.iterations = iterations;
        this.bodies = [];
        this.contacts = [];
    }

    step() {
        this.contacts.length = 0;
        for (let i = 0; i < this.bodies.length; ++i) {
            let bodyA = bodies[i];
            for(let j = i + 1; j < this.bodies.length; ++j) {
                let bodyB = bodies[j];

                if(bodyA.invMass === 0 && bodyB.invMass === 0) {
                    continue;
                }

                /*Manifold m = new Manifold(bodyA, bodyB);
                m.solve();

                if (m.contactCount > 0) {
                    this.contacts.add(m);
                }*/
            }
        }

        /*
        // Integrate forces
        for (let i = 0; i < this.bodies.length; ++i) {
            this.integrateForces(this.bodies[i], this.dt);
        }

        // Initialize collision
        for (let i = 0; i < this.bodies.length; ++i) {
            this.contacts[i].initialize();
        }

        // Solve collisions
        for (let j = 0; j < this.iterations; ++j) {
            for (let i = 0; i < this.contacts.length; ++i) {
                this.contacts[i].applyImpulse();
            }
        }

        // Integrate velocities
        for (let i = 0; i < this.bodies.length; ++i) {
            this.integrateVelocity(this.bodies[i], dt);
        }

        // Correct positions
        for (let i = 0; i < this.contacts.length; ++i) {
            this.contacts[i].positionalCorrection();
        }

        // Clear all forces
        for (let i = 0; i < this.bodies.length; ++i) {
            let b = bodies[i];
            b.force.set(0, 0);
            b.torque = 0;
        }
        */
    }

    add(shape, x, y) {
        let b = new Body(shape, x, y);
        this.bodies.push(b);
        return b;
    }

    clear() {
        this.contacts.length = 0;
        this.bodies.length = 0;
    }

    integrateForces(body, dt) {
        if (body.invMass === 0) {
            return;
        }

        let dts = dt * 0.5;
        body.velocity.addMult(body.force, body.invMass * dts);
        body.velocity.addMult(ImpulseMath.GRAVITY, dts);
        body.angularVelocity += body.torque * body.invInertia * dts;
    }

    integrateVelocity(body, dt) {
        if(body.invMass === 0) {
            return;
        }

        body.position.addMult(body.velocity, dt);
        body.orient += body.angularVelocity * dt;
        body.setOrient(body.orient);

        this.integrateForces(body, dt);
    }

    render(canvas, ctx) {
        ctx.clearRect(0, 0, canvas.width, canvas.height);

        // Render bodies
        for (let i = 0; i < this.bodies.length; ++i) {
            this.bodies[i].renderBody(ctx);
        }
    }
}

let scene = new ImpulseScene(0, 0);
scene.add(new Shape({
    type: 'c',
    radius: 10
}), 100, 100);
scene.add(new Shape({
    type: 'c',
    radius: 100
}), 500, 300);

function loop() {
    scene.render(canvas, ctx);
    requestAnimationFrame(loop);
}
loop();