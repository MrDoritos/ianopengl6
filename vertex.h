#pragma once
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <GL/glew.h>

#define BOTTOM_FACE 4
#define RIGHT_FACE 0
#define LEFT_FACE 2
#define BACK_FACE 3
#define FRONT_FACE 1
#define TOP_FACE 5

typedef unsigned int guint;
typedef unsigned char guint8;

struct vec3d : public glm::vec3 {
	vec3d(float x, float y, float z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}
	vec3d() {
		this->x = 0;
		this->y = 0;
		this->z = 0;
	}
    vec3d operator+(vec3d b) {
        return vec3d(x + b.x, y + b.y, z + b.z);
    }

    vec3d operator-(vec3d b) {
        return vec3d(x - b.x, y - b.y, z - b.z);
    }
	vec3d(float xyz):vec3d(xyz,xyz,xyz){}
};

GLfloat textureCube3[] = {	
    // Front face
    0.0f, 1.0f, // top-right
    1.0f, 0.0f, // bottom-right      
    0.0f, 0.0f, // bottom-left  
	
    0.0f, 0.0f, // bottom-left
    1.0f, 0.0f, // top-left       
    1.0f, 1.0f, // top-right
	
    // Right face
    1.0f, 0.0f, // top-left
    1.0f, 1.0f, // top-right      
    0.0f, 1.0f, // bottom-right          
    0.0f, 1.0f, // bottom-right
    0.0f, 0.0f, // bottom-left
    1.0f, 0.0f, // top-left
	
    // Back face
    0.0f, 0.0f, // Bottom-left
    1.0f, 0.0f, // bottom-right    
    1.0f, 1.0f, // top-right              
    1.0f, 1.0f, // top-right
    0.0f, 1.0f, // top-left
    0.0f, 0.0f, // bottom-left                
	
    // Left face
    0.0f, 1.0f, // bottom-left
    1.0f, 0.0f, // top-right
    1.0f, 1.0f, // top-left       
	
    0.0f, 1.0f, // bottom-left
    0.0f, 0.0f, // bottom-right
    1.0f, 0.0f, // top-right
	
	
    // Bottom face          
    0.0f, 1.0f, // top-right
    1.0f, 1.0f, // top-left        
    1.0f, 0.0f, // bottom-left
	
    1.0f, 0.0f, // bottom-left
    0.0f, 0.0f, // bottom-right
    0.0f, 1.0f, // top-right
	
    // Top face
    0.0f, 1.0f, // top-left
    1.0f, 1.0f, // top-right
    1.0f, 0.0f, // bottom-right                 
    1.0f, 0.0f, // bottom-right
    0.0f, 0.0f, // bottom-left  
    0.0f, 1.0f  // top-left              	
};

GLfloat positionCube3[] = {	
	//Front BAD
	1.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f,
	0.0f, 1.0f, 0.0f,
	
	1.0f, 0.0f, 0.0f,
	0.0f, 1.0f, 0.0f,
	1.0f, 1.0f, 0.0f,
	
	//Right GOOD
	1.0f, 0.0f, 1.0f,
	1.0f, 0.0f, 0.0f,
	1.0f, 1.0f, 0.0f,
	
	1.0f, 0.0f, 1.0f,
	1.0f, 1.0f, 0.0f,
	1.0f, 1.0f, 1.0f,
	
	//Back GOOD
	0.0f, 0.0f, 1.0f,
	1.0f, 0.0f, 1.0f,
	1.0f, 1.0f, 1.0f,
	
	0.0f, 0.0f, 1.0f,
	1.0f, 1.0f, 1.0f,
	0.0f, 1.0f, 1.0f,
			
	//Left GOOD
	0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 1.0f,
	0.0f, 1.0f, 1.0f,
	
	0.0f, 0.0f, 0.0f,
	0.0f, 1.0f, 1.0f,
	0.0f, 1.0f, 0.0f,
	
	//Bottom BAD
	1.0f, 0.0f, 0.0f,
	1.0f, 0.0f, 1.0f,
	0.0f, 0.0f, 1.0f,
	
	1.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 1.0f,
	0.0f, 0.0f, 0.0f,
	
	//Top GOOD
	0.0f, 1.0f, 1.0f,
	1.0f, 1.0f, 1.0f,
	1.0f, 1.0f, 0.0f,
	
	0.0f, 1.0f, 1.0f,
	1.0f, 1.0f, 0.0f,
	0.0f, 1.0f, 0.0f,
};

int map[2][3][2] = {
    //New order 0 - 1 - 2  2 - 3 - 0
    //Old order 0 - 2 - 1  3 - 1 - 0
    { //Triangle 0
        { //Vertex 0
            0, 3
        }, 
        { //Vertex 1
            2, 3
        }, 
        { //Vertex 2
            2, 1
        }
    },
    { //Triangle 1
        { //Vertex 0
            0, 3
        },
        { //Vertex 1
            2, 1
        },
        { //Vertex 2
            0, 1
        }
    }			
};