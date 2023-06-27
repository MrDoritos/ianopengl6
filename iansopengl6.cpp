#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <fstream>
#include <sstream>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

typedef unsigned int guint;
typedef unsigned char guint8;

struct world_t;
struct game_t;
struct blockaccessor_t;
struct block_t;
struct fullblock_t;
struct blockstate_t;
struct chunk_t;

glm::mat4 proj;
glm::mat4 view;
float startYOffset = 0.0f;
glm::vec3 cameraPos = glm::vec3(0.0f, 0.0f + startYOffset + 0.0f, 0.0f);
glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
glm::vec3 cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);
bool firstMouse = true;
float yaw = -90.0f;
float pitch = 0.0f;
float lastX = 800.0f / 2.0;
float lastY = 600.0 / 2.0;
bool showTriangles = false;
bool doDebug = true;
float sunAngle = 0.5f;
int ticks = 0;
int selectedBlock = 0;
float deltaTime = 0.0f;
float lastFrame = 0.0f;
GLint uniModel;
GLint uniView;
GLuint shaderProgram;
GLFWmonitor* monitor = glfwGetPrimaryMonitor();
void key_callback(GLFWwindow* window);
void key_frame(GLFWwindow* window);
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);

class perlin {
	public:
		static float getPerlin(float x, float y);
		static float getPerlin(float x, float y, int octaves, float persistence);
		static float getNormal(float omax, float omin, float max, float min, float value);
		static float getNormalNoise(float x, float z);
		static float cosine_interpolate(float a, float b, float x);
		static float smooth_noise_2D(float x, float y);
		static float interpolated_noise(float x, float y);
		static float noise(int x, int y);
		
		static float octaves;
		static float persistence;
		static float lacunarity;
};

#define PI 3.1415927

float perlin::noise(int x, int y) {
    int n = x + y * 57;
    n = (n<<13) ^ n;
    return (1.0 - ( (n * ((n * n * 15731) + 789221) +  1376312589) & 0x7fffffff) / 1073741824.0);
}

float perlin::cosine_interpolate(float a, float b, float x) {
    float ft = x * PI;
    float f = (1 - cos(ft)) * 0.5;
    float result =  a*(1-f) + b*f;
    return result;
}

float perlin::smooth_noise_2D(float x, float y) {  
    float corners = ( noise(x-1, y-1)+noise(x+1, y-1)+noise(x-1, y+1)+noise(x+1, y+1) ) / 16;
    float sides   = ( noise(x-1, y)  +noise(x+1, y)  +noise(x, y-1)  +noise(x, y+1) ) /  8;
    float center  =  noise(x, y) / 4;

    return corners + sides + center;
}

float perlin::interpolated_noise(float x, float y) {
    int x_whole = (int) x;
    float x_frac = x - x_whole;

    int y_whole = (int) y;
    float y_frac = y - y_whole;

    float v1 = smooth_noise_2D(x_whole, y_whole); 
    float v2 = smooth_noise_2D(x_whole, y_whole+1); 
    float v3 = smooth_noise_2D(x_whole+1, y_whole); 
    float v4 = smooth_noise_2D(x_whole+1, y_whole+1); 

    float i1 = cosine_interpolate(v1,v3,x_frac);
    float i2 = cosine_interpolate(v2,v4,x_frac);

    return cosine_interpolate(i1, i2, y_frac);
}

float perlin::octaves = 8.0f;
float perlin::persistence = 0.5f;
float perlin::lacunarity = 1.0f;

float perlin::getPerlin(float x, float y, int octaves, float persistence) {
    float total = 0;

    for(int i=0; i<octaves-1; i++)
    {
        float frequency = pow(2,i);
        float amplitude = pow(persistence,i);
        total = total + interpolated_noise(x * frequency, y * frequency) * amplitude;
    }
    return total;
}

float perlin::getPerlin(float x, float y) {
    float total = 0;
	//x += int(1 << 20);
	//y += int(1 << 20);

    for(int i=0; i<octaves-1; i++)
    {
        float frequency = pow(2,i);
        float amplitude = pow(persistence,i);
        total = total + interpolated_noise(x * frequency, y * frequency) * amplitude;
    }
    return total;
}
 
float perlin::getNormal(float omax, float omin, float max, float min, float value) {
	return (max - min) / (omax - omin) * (value - omax) + max;
}

float perlin::getNormalNoise(float x, float z) {
	float noise = perlin::getPerlin(x,z);
	return perlin::getNormal(1,-1,1,0,noise);
}

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
    // Back face
    0.0f, 0.0f, 0.0f,
    1.0f, 0.0f, 0.0f,
    1.0f, 1.0f, 0.0f,
    1.0f, 1.0f, 0.0f,
    0.0f, 1.0f, 0.0f,
    0.0f, 0.0f, 0.0f,
    // Front face
    0.0f, 0.0f, 1.0f,
    1.0f, 1.0f, 1.0f,
    1.0f, 0.0f, 1.0f,
    1.0f, 1.0f, 1.0f,
    0.0f, 0.0f, 1.0f,
    0.0f, 1.0f, 1.0f,
    // Left face
    0.0f, 1.0f, 1.0f,
    0.0f, 0.0f, 0.0f,
    0.0f, 1.0f, 0.0f,
    0.0f, 0.0f, 0.0f,
    0.0f, 1.0f, 1.0f,
    0.0f, 0.0f, 1.0f,
    // Right face
    1.0f, 1.0f, 1.0f,
    1.0f, 1.0f, 0.0f,
    1.0f, 0.0f, 0.0f,
    1.0f, 0.0f, 0.0f,
    1.0f, 0.0f, 1.0f,
    1.0f, 1.0f, 1.0f,
    // Bottom face          
    0.0f, 0.0f, 0.0f,
    1.0f, 0.0f, 1.0f,
    1.0f, 0.0f, 0.0f,
    1.0f, 0.0f, 1.0f,
    0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 1.0f,
    // Top face
    0.0f, 1.0f, 0.0f,
    1.0f, 1.0f, 0.0f,
    1.0f, 1.0f, 1.0f,
    1.0f, 1.0f, 1.0f,
    0.0f, 1.0f, 1.0f,
    0.0f, 1.0f, 0.0f
};

void errorCallback(int error, const char* description) {
    std::cout << "GLFW Error: " << description << std::endl;
}

struct game_t {
    struct world_t *currentWorld;
    std::vector<vec3d> vectorData;
    void render();
    void mesh(chunk_t *chunk);
};

struct blockstate_t {
    short blockId;
};

struct block_t {
    static std::vector<block_t*> blocks;

    static block_t *air;
    static block_t *stone;

    block_t() {
        blockId = 0;
        
    }

    block_t(int blockId) {
        this->blockId = blockId;
        blocks.push_back(this);
    }


    virtual blockstate_t getDefaultState() {
        blockstate_t ret;
        ret.blockId = blockId;
        return ret;
    }

    int blockId;  
};

std::vector<block_t*> block_t::blocks;
block_t *block_t::air;
block_t *block_t::stone;

struct air : public block_t {
    air():block_t(0){}
};

struct stone : public block_t {
    stone():block_t(1){}
};

struct fullblock_t {
    fullblock_t() {
        state = nullptr;
        block = nullptr;
    }
    vec3d position;
    blockstate_t *state;
    block_t *block;
};

struct blockaccessor_t {
    virtual fullblock_t getBlockAtAbsolute(vec3d pos) = 0;
    virtual fullblock_t getBlockAtRelative(vec3d pos) = 0;
    virtual bool isAvailable() = 0;
    virtual bool isBlockLoaded(vec3d pos) = 0;
};

static union _draw_type {
    struct { float x,y,z,u,v,i,a; } a;
    float p[7];
} *draw_buffer = 0;
                    
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
            2, 1
        },
        { //Vertex 1
            0, 1
        },
        { //Vertex 2
            0, 3
        }
    }			
};

game_t currentGame;
int drawen_buffers;

struct octree_chunk_t : public blockaccessor_t {
    fullblock_t getBlockAtAbsolute(vec3d pos) override {
        int width = std::cbrt(pow(8, size));
        int halfwidth = width;

        vec3d pos_rel = pos - position;

        fullblock_t nothing_to_return;
        static blockstate_t *nothing_blockstate = new blockstate_t();
        nothing_blockstate->blockId = -1;
        nothing_to_return.block = block_t::air;
        nothing_to_return.state = nothing_blockstate;

        if (pos.x < position.x ||
            pos.x >= position.x + halfwidth ||
            pos.y < position.y||
            pos.y >= position.y + halfwidth ||
            pos.z < position.z ||
            pos.z >= position.z + halfwidth) {
                return nothing_to_return;
            }

        if (size > 2) {
            int childx = (pos_rel.x)/width;
            int childy = (pos_rel.y)/width;
            int childz = (pos_rel.z)/width;
            return children[(childz * 4) + (childy * 2) + childx].getBlockAtAbsolute(pos);
        } else {
            if (!point) {
                //puts("Nothing here");
                return nothing_to_return;
            }

            int x = int(pos_rel.x);
            int y = int(pos_rel.y);
            int z = int(pos_rel.z);
            blockstate_t *state = &point[(y * 16) + (z * 4) + x];
            fullblock_t ret;
            ret.state = state;
            if (state->blockId < block_t::blocks.size())
                ret.block = block_t::blocks.at(state->blockId);
            else
                ret.block = block_t::air;
            ret.position = pos;
            return ret;
        }

        return nothing_to_return;
    }

    fullblock_t getBlockAtRelative(vec3d pos) override {
        return getBlockAtAbsolute(pos + position);
    }

    bool isAvailable() override {
        return true;
    }

    bool isBlockLoaded(vec3d pos) override {
        return true;
    }

    octree_chunk_t() {
        has_children = false;
        is_meshable = false;
        needsMeshed = true;
        children = 0;
        point = 0;
        vbo = 0;
        vao = 0;
        info_count = 0;
        this->position = vec3d{0,0,0};
    }

    octree_chunk_t(int _size, vec3d _position) : size(_size) {
        has_children = false;
        is_meshable = false;
        needsMeshed = true;
        children = 0;
        point = 0;
        vbo = 0;
        vao = 0;
        info_count = 0;
        this->position = _position;

        if (size >= 2) { //If 4x4x4 and above, we hold blocks and can render individually
            glGenBuffers(1, &vbo);
            glGenVertexArrays(1, &vao);
            is_meshable = true;
        }
        if (size == 2) { //If 4x4x4, we hold blocks
            point = new blockstate_t[64];
            memset(point, 0, sizeof(*point) * 64);
            blockstate_t state = block_t::stone->getDefaultState();
            //point[0] = block_t::stone->getDefaultState();
            //point[0] = block_t::stone->getDefaultState();
            //point[0] = block_t::stone->getDefaultState();
            //point[0] = block_t::stone->getDefaultState();
            //point[0] = state;
            //point[3] = state;
            //point[12] = state;
            //point[15] = state;
            //point[48] = state;
            //point[51] = state;
            //point[60] = state;
            //point[63] = state;
            generate();
        }
        if (size > 2) { //All sizes above 4x4x4 create children
            generate_children();
        }
    }

    void generate() {
        for (int y = 0; y < 4; y++)
            for (int x = 0; x < 4; x++)
                for (int z = 0; z < 4; z++) {
                    if (position.y + y > 0 && position.y + y < 5)
                    point[(y * 16) + (z * 4) + x] = block_t::stone->getDefaultState();
                }
                    //*getBlockAtRelative(vec3d{x,y,z}).state = block_t::stone->getDefaultState();
//return;
        for (int x = 0; x < 4; x++) {
            for (int z = 0; z < 4; z++) {
                int ofx = x + position.x;
                int ofz = z + position.z;
                perlin::octaves = 2;
                perlin::persistence = 0.55f;
                //float scale = 0.005f;
                float scale = 0.01f;
                float height = 80.0f;
                float value = perlin::getPerlin(ofx * scale + 16000.0f, ofz * scale + 16000.0f) * height;

                for (int y = 0; y < 4; y++) {
                    if (value - 1.0f > position.y + y)  
                        point[(y * 16) + (z * 4) + x] = block_t::stone->getDefaultState();
                }
            }
        }
    }

    GLuint vbo, vao; //size >= 5 32x32x32 or more
    blockstate_t *point; //size = 3 4x4x4
    const int size;
    octree_chunk_t *children; //size > 1
    bool has_children;
    bool is_meshable;
    bool needsMeshed;
    int info_count;
    vec3d position;
    static const int info_max = 800000;

    bool needsToRender() {
        if (!is_meshable || !children || size < 3)
            return false;
        //return true;
        float scale = 0.1f;
        float threshold = scale * size;
        threshold = std::cbrt(pow(8, size)) * 4;
        float dist = 
        sqrtf(((cameraPos.x - position.x) * (cameraPos.x - position.x)) +
              ((cameraPos.y - position.y) * (cameraPos.y - position.y)) +
              ((cameraPos.z - position.z) * (cameraPos.z - position.z)));

        //return (dist < pow(2, size-1));
        return dist < threshold;
    }

    void generate_children() {
        printf("Generating size %i\n", size);
        {
            double childSize = std::cbrt(pow(8, size-1));
            double halfSize = childSize;//4.0f;// / 2.0;
            
            children = (octree_chunk_t*)malloc(sizeof(octree_chunk_t)*8);
            //printf("halfsize %f\n", halfSize);
            // Generate the eight children origins based on the position and size of the parent node
            /*
            new (&children[0]) octree_chunk_t(size-1, vec3d{ position.x - halfSize, position.y - halfSize, position.z - halfSize });
            new (&children[1]) octree_chunk_t(size-1, vec3d{ position.x + halfSize, position.y - halfSize, position.z - halfSize });
            new (&children[2]) octree_chunk_t(size-1, vec3d{ position.x - halfSize, position.y + halfSize, position.z - halfSize });
            new (&children[3]) octree_chunk_t(size-1, vec3d{ position.x + halfSize, position.y + halfSize, position.z - halfSize });
            new (&children[4]) octree_chunk_t(size-1, vec3d{ position.x - halfSize, position.y - halfSize, position.z + halfSize });
            new (&children[5]) octree_chunk_t(size-1, vec3d{ position.x + halfSize, position.y - halfSize, position.z + halfSize });
            new (&children[6]) octree_chunk_t(size-1, vec3d{ position.x - halfSize, position.y + halfSize, position.z + halfSize });
            new (&children[7]) octree_chunk_t(size-1, vec3d{ position.x + halfSize, position.y + halfSize, position.z + halfSize });
            */
            new (&children[0]) octree_chunk_t(size-1, vec3d{ position.x, position.y, position.z });
            new (&children[1]) octree_chunk_t(size-1, vec3d{ position.x + halfSize, position.y, position.z });
            new (&children[2]) octree_chunk_t(size-1, vec3d{ position.x, position.y + halfSize, position.z });
            new (&children[3]) octree_chunk_t(size-1, vec3d{ position.x + halfSize, position.y + halfSize, position.z });
            new (&children[4]) octree_chunk_t(size-1, vec3d{ position.x, position.y, position.z + halfSize });
            new (&children[5]) octree_chunk_t(size-1, vec3d{ position.x + halfSize, position.y, position.z + halfSize });
            new (&children[6]) octree_chunk_t(size-1, vec3d{ position.x, position.y + halfSize, position.z + halfSize });
            new (&children[7]) octree_chunk_t(size-1, vec3d{ position.x + halfSize, position.y + halfSize, position.z + halfSize });

            has_children = true;
        }
    }

    ~octree_chunk_t() {
        glDeleteBuffers(1, &vbo);
        glDeleteVertexArrays(1, &vao);
    }

    void mesh(_draw_type *buffer, int *index);

    void render() {
        //if (!is_meshable)
        //    return;

        //if (!chunk->loaded) {
        //    fprintf(stderr, "Not yet loaded\n");
        //    continue;
        //}

        bool childRendered = false;
        for (int i = 0; i < 8; i++) {
            if (children[i].needsToRender()) {
                children[i].render();
                childRendered = true;
            }
        }
        if (childRendered) //For some reason this turns on/off lod
            return;//return;//return; //If a child rendered, a higher LOD is available and was used, making this size/LOD useless


        if (needsMeshed) {
            if (!draw_buffer) {
                printf("Allocating %i entries\n", info_max);
                draw_buffer = (_draw_type*)malloc(sizeof(_draw_type) * info_max);
            }
            needsMeshed = false;
            info_count = 0;

            fprintf(stderr, "Meshing\n");            
            mesh(draw_buffer, &info_count);

            if (info_count == 0)
                return;

            glBindVertexArray(vao);
            glBindBuffer(GL_ARRAY_BUFFER, vbo);

            glBufferData(GL_ARRAY_BUFFER, info_count * sizeof(*draw_buffer), draw_buffer, GL_STATIC_DRAW);

            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 28, 0);
            glEnableVertexAttribArray(0);

            glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 28, 12);
            glEnableVertexAttribArray(1);

            glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 28, 20);
            glEnableVertexAttribArray(2);

            needsMeshed = false;
            printf("Uploaded %i verticies\n", info_count);

        }


        if (info_count == 0)
            return;
        //glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glBindVertexArray(vao);
        //glBindTexture(GL_TEXTURE_2D, textureId);

        glm::mat4 model = glm::translate(glm::mat4(1.0), vec3d{0,0,0});
        glUniformMatrix4fv(uniModel, 1, GL_FALSE, glm::value_ptr(model));

        // Draw the quad
        //glBegin(GL_LINES);
        //printf("Drawing %i verticies\n", chunk->info_count);
        glDrawArrays(GL_TRIANGLES, 0, info_count);
        drawen_buffers++;
        //glEnd();
        // Unbind the VAO
        glBindVertexArray(0);

    }
};

struct chunk_t : public blockaccessor_t {
    bool loaded, needsMeshed;
    GLuint vbo, vao;
    int info_count;
    vec3d position;
    blockstate_t chunkData[16][256][16];

    chunk_t() {
        loaded = false;
        needsMeshed = false;
    }

    fullblock_t getBlockAtAbsolute(vec3d pos) override {
        return getBlockAtRelative(vec3d{int(pos.x)%16,int(pos.y),int(pos.z)%16});
    }
    
    fullblock_t getBlockAtRelative(vec3d pos) override {
        fullblock_t ret;
        ret.state = &chunkData[int(pos.x)][int(pos.y)][int(pos.z)];
        ret.position = pos + position;
        ret.block = block_t::blocks.at(ret.state->blockId);
        return ret;
    }

    bool isAvailable() override {
        return loaded;
    }

    bool isBlockLoaded(vec3d pos) override {
        return (isAvailable() &&
                pos.x >= position.x && pos.x < position.x + 16 &&
                pos.y >= position.y && pos.y < position.y + 256 &&
                pos.z >= position.z && pos.z < position.z + 16);
    }

    void generate(vec3d newPosition) {
        puts("Generating");
        memset(&chunkData[0], 0, sizeof(blockstate_t) * 16 * 256 * 16);
        for (int x = 0; x < 16; x++)
        for (int z = 0; z < 16; z++)
            chunkData[x][16][z] = block_t::blocks.at(1)->getDefaultState();

        position = newPosition;
        loaded = true;
        needsMeshed = true;
    }
};

struct world_t : public blockaccessor_t {
    //std::vector<chunk_t*> chunks;
    octree_chunk_t *chunks = new octree_chunk_t(8, vec3d{0,0,0});

    fullblock_t getBlockAtAbsolute(vec3d pos) override {
        return chunks->getBlockAtAbsolute(pos);
    }

    fullblock_t getBlockAtRelative(vec3d pos) override {
        //for (chunk_t* chunk : chunks) {
        //    if (chunk->isBlockLoaded(pos))
        //        return chunk->getBlockAtAbsolute(pos);
        //}

        return chunks->getBlockAtRelative(pos);
    }

    bool isBlockLoaded(vec3d pos) override {
        //for (auto chunk : chunks)
        //    if (chunk->isBlockLoaded(pos))
        //        return true;
        //return false;
        return chunks->isBlockLoaded(pos);
    }

    bool isAvailable() override {
        return true;
    }
};

void octree_chunk_t::mesh(_draw_type* buffer, int* index) {
        int ofx = position.x;
        int ofy = position.y;
        int ofz = position.z;
        //Wireframe around chunk


                        for (int l = 0; l < 12; l++) {
                            for (int l1 = 0; l1 < 3; l1++) {			
                                if ((*index) + 2 > info_max) {
                                    puts("Could not resize buffers");
                                    //exit(1);b
                                    printf("%i/%i %p\n", *index, info_max, index);
                                    break;
                                }
                                break;
                                float factor = (1.0f/16.0f);
                                float textureAtlas[] = {0*factor,0*factor,1*factor,1*factor};
                                //SHOWS ONLY ONES INDIVIDUALLY RENDERED BECAUSE OF APOSITION SHADER
                                draw_buffer[*index].a.x = positionCube3[(l * 9) + (l1 * 3) + 0] + ofx;// + (float(ofx) / 128.0f);
                                draw_buffer[*index].a.y = positionCube3[(l * 9) + (l1 * 3) + 1] + ofy;// + (float(ofy) / 128.0f);
                                draw_buffer[*index].a.z = positionCube3[(l * 9) + (l1 * 3) + 2] + ofz;// + (float(ofz) / 128.0f);
                                
                                draw_buffer[*index].a.u = textureAtlas[map[l % 2][l1][0]];
                                draw_buffer[*index].a.v = textureAtlas[map[l % 2][l1][1]];
                                
                                
                                draw_buffer[*index].a.i = 1.0f;
                                draw_buffer[*index].a.a = 1.0f;
                                
                                (*index)++;
                            }
                        }

        //If we have children that may use their vbos, we either add them to our mesh or let them handle it
        if (size > 2) {
            bool childNotMeshed = true;
            if (has_children) {
                for (int i = 0; i < 8; i++) {
                    octree_chunk_t* child = &children[i];
                    if (!child->needsToRender()) {
                        child->mesh(buffer, index); //Nothing renders if this is commented out
                    } else {
                        childNotMeshed = false;
                    }
                }
            }

            if (!childNotMeshed)
                return; //If no child contributed to the mesh, we need to draw ourselves
        }

        if (size > 2)
            return;

        //We have blocks that can be meshed
        //puts("Mesh");
        //printf("%i %i %i\n", ofx, ofy, ofz);

        int width = std::cbrt(pow(8, size));
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < width; y++) {
                for (int z = 0; z < width; z++) {
                    
                    //getBlockAt(x + ofx, y + ofy, z + ofz, &block);
                    //if (block.parent == blocks::AIR)
                    //	continue;
                    //printf("[%i,%i,%i]:%i ", x,y,z,block.parent->id);
                    fullblock_t fb = getBlockAtRelative(vec3d{x,y,z});
                    if (fb.block == block_t::air)
                        continue;

                    auto _outOfChunkCheck = [&](int _x, int _y, int _z) {
                        //return rand() % 3 == 0;
                        //return true;
                        //return false;
                        //puts("Out of chunk check");
                        fullblock_t _fb = currentGame.currentWorld->getBlockAtAbsolute(vec3d{_x,_y,_z});
                        if (_fb.state == nullptr)
                            return true;
                        return _fb.state->blockId == 0;
                    };
                    
                    auto isVisible = [&](int _x, int _y, int _z) {
                        //_x += ofx; _y += ofy; _z += ofz;
                        
                        if (_x > width-1 ||
                            _y > width-1 ||
                            _z > width-1)
                            return _outOfChunkCheck(_x + ofx, _y + ofy, _z + ofz);

                        if (_x < 0 ||
                            _y < 0 ||
                            _z < 0)
                            return _outOfChunkCheck(_x + ofx, _y + ofy, _z + ofz);
                        
                        //return true;
                        fullblock_t blk = getBlockAtRelative(vec3d{_x, _y, _z});
                        return blk.state->blockId == 0;
                    };

                    auto _addFace = [&](int _ofs) {
                        
                        for (int l = _ofs; l < _ofs + 2; l++) {
                            for (int l1 = 0; l1 < 3; l1++) {			
                                if ((*index) + 2 > info_max) {
                                    puts("Could not resize buffers");
                                    //exit(1);b
                                    printf("%i/%i %p\n", *index, info_max, index);
                                    break;
                                }
                            
                                float factor = (1.0f/16.0f);
                                float textureAtlas[] = {7*factor,0*factor,8*factor,1*factor};

                                draw_buffer[*index].a.x = positionCube3[(l * 9) + (l1 * 3) + 0] + float(x + ofx);
                                draw_buffer[*index].a.y = positionCube3[(l * 9) + (l1 * 3) + 1] + float(y + ofy);
                                draw_buffer[*index].a.z = positionCube3[(l * 9) + (l1 * 3) + 2] + float(z + ofz);
                                
                                draw_buffer[*index].a.u = textureAtlas[map[l % 2][l1][0]];
                                draw_buffer[*index].a.v = textureAtlas[map[l % 2][l1][1]];
                                
                                
                                draw_buffer[*index].a.i = 1.0f;
                                draw_buffer[*index].a.a = 1.0f;
                                
                                (*index)++;
                            }
                        }
                        
                    };
                    
                    
                    //Old order  10 -> 8 -> 0 -> 4 -> 6 -> 2
                    //New order? 10 -> 8 -> 6 -> 4 -> 0 -> 2
                            if (isVisible(x, y - 1, z)) //Bottom 2
                                _addFace(8);  //FINE
                            if (isVisible(x, y, z - 1)) //Right 4
                                _addFace(0);//_addFace(0);  //FINE					
                            if (isVisible(x, y, z + 1)) //Left 6
                                _addFace(4);//_addFace(4);  //FINE				
                            if (isVisible(x - 1, y, z)) //Back 8
                                _addFace(6);//_addFace(6); //FINE						
                            if (isVisible(x + 1, y, z)) //Front 10
                                _addFace(2); //FINE
                            if (isVisible(x, y + 1, z))
                                _addFace(10);

                }
            }
        }



        //size == 2, mesh blocks
}


void start_game() {
    puts("Starting game");

    block_t::air = new air;
    block_t::stone = new stone;

    currentGame.currentWorld = new world_t();
    /*
    for (int x = 0; x < 3; x++) {
        for (int z = 0; z < 3; z++) {
            chunk_t *newChunk = new chunk_t();
            newChunk->generate(vec3d{x * 16,0,z * 16});
            glGenBuffers(1, &newChunk->vbo);
            glGenVertexArrays(1, &newChunk->vao);
            currentGame.currentWorld->chunks.push_back(newChunk);
        }
    }
    */
}

void game_t::mesh(chunk_t* chunk) {
    glBindVertexArray(chunk->vao);
    glBindBuffer(GL_ARRAY_BUFFER, chunk->vbo);
    float quadpoints[] = {
     //0.0f,  0.0f,  0.0f,     0.5f, 0.5f,    1.0f, 1.0f,
     0.0f, 0.0f,  0.0f,     0.0f, 0.0f,    1.0f, 1.0f,
     1.0f, 0.0f,  0.0f,     1.0f, 0.0f,    1.0f, 1.0f,
     0.0f,  1.0f,  0.0f,     0.0f, 1.0f,    1.0f, 1.0f,
     1.0f,  1.0f,  0.0f,     1.0f, 1.0f,    1.0f, 1.0f
    };

    static float quad_front[] = {
        0.0f, 0.0f, 1.0f,
        1.0f, 1.0f, 1.0f,
        1.0f, 0.0f, 1.0f,
        1.0f, 1.0f, 1.0f,
        0.0f, 0.0f, 1.0f,
        0.0f, 1.0f, 1.0f
    };

    static float text_front[] = {
        1.0f, 0.0f, // top-left
        1.0f, 1.0f, // top-right      
        0.0f, 1.0f, // bottom-right          
        0.0f, 1.0f, // bottom-right
        0.0f, 0.0f, // bottom-left
        1.0f, 0.0f // top-left
    };

    static const GLfloat cube_strip[] = {
    0.f, 1.f, 1.f,     // Front-top-left
    1.f, 1.f, 1.f,      // Front-top-right
    0.f, 0.f, 1.f,    // Front-bottom-left
    1.f, 0.f, 1.f,     // Front-bottom-right
    1.f, 0.f, 0.f,    // Back-bottom-right
    1.f, 1.f, 1.f,      // Front-top-right
    1.f, 1.f, 0.f,     // Back-top-right
    0.f, 1.f, 1.f,     // Front-top-left
    0.f, 1.f, 0.f,    // Back-top-left
    0.f, 0.f, 1.f,    // Front-bottom-left
    0.f, 0.f, 0.f,   // Back-bottom-left
    1.f, 0.f, 0.f,    // Back-bottom-right
    0.f, 1.f, 0.f,    // Back-top-left
    1.f, 1.f, 0.f      // Back-top-right
    };

    int info_max = 10000;

    if (!draw_buffer) {
        printf("Allocating %i entries\n", info_max);
        draw_buffer = (_draw_type*)malloc(sizeof(_draw_type) * info_max);
    }

    {
        //float draw_info[info_max]; // 400 faces for testing
        chunk->info_count = 0;

        puts("Mesh");
        int ofx = chunk->position.x * 16;
        int ofy = chunk->position.y * 16;
        int ofz = chunk->position.z * 16;

        for (int x = 0; x < 16; x++) {
            for (int y = 0; y < 256; y++) {
                for (int z = 0; z < 16; z++) {
                    
                    //getBlockAt(x + ofx, y + ofy, z + ofz, &block);
                    //if (block.parent == blocks::AIR)
                    //	continue;
                    //printf("[%i,%i,%i]:%i ", x,y,z,block.parent->id);
                    fullblock_t fb = chunk->getBlockAtRelative(vec3d{x,y,z});
                    if (fb.block == block_t::air)
                        continue;
                    auto _outOfChunkCheck = [&](int _x, int _y, int _z) {
                        //return true;
                        fullblock_t _fb = currentGame.currentWorld->getBlockAtAbsolute(vec3d{_x,_y,_z});
                        if (_fb.state == nullptr)
                            return false;
                        return _fb.block == block_t::air;
                    };
                    
                    auto isVisible = [&](int _x, int _y, int _z) {
                        if (_x < 0 || _y < 0 || _z < 0)
                            return _outOfChunkCheck(_x + x, _y + y, _z + z);
                        if (_x > 15 ||
                            _y > 255 ||
                            _z > 15)
                            return _outOfChunkCheck(_x + x, _y + y, _z + z);
                        
                        
                        fullblock_t blk = chunk->getBlockAtAbsolute(vec3d{_x + ofx, _y + ofy, _z + ofz});
                        return blk.block == block_t::air;
                    };
                    
                    auto _getBlockId = [&](int _x, int _y, int _z) {
                        return chunk->getBlockAtAbsolute(vec3d{_x + ofx, _y + ofy, z + ofz}).block->blockId;
                    };

                    auto _addFace = [&](int _ofs) {
                        
                        for (int l = _ofs; l < _ofs + 2; l++) {
                            for (int l1 = 0; l1 < 3; l1++) {			
                                if (chunk->info_count + 2 > info_max) {
                                    puts("Could not resize buffers");
                                    //exit(1);b
                                    printf("%i/%i\n", chunk->info_count, info_max);
                                    break;
                                }
                            
                                float factor = (1.0f/16.0f);
                                float textureAtlas[] = {1*factor,0*factor,2*factor,1*factor};

                                draw_buffer[chunk->info_count].a.x = positionCube3[(l * 9) + (l1 * 3) + 0] + float(x);
                                draw_buffer[chunk->info_count].a.y = positionCube3[(l * 9) + (l1 * 3) + 1] + float(y);
                                draw_buffer[chunk->info_count].a.z = positionCube3[(l * 9) + (l1 * 3) + 2] + float(z);
                                
                                draw_buffer[chunk->info_count].a.u = textureAtlas[map[l % 2][l1][0]];
                                draw_buffer[chunk->info_count].a.v = textureAtlas[map[l % 2][l1][1]];
                                
                                
                                draw_buffer[chunk->info_count].a.i = 1.0f;
                                draw_buffer[chunk->info_count].a.a = 1.0f;
                                
                                chunk->info_count++;
                            }
                        }
                        
                    };
                    
                    
                    //Old order  10 -> 8 -> 0 -> 4 -> 6 -> 2
                    //New order? 10 -> 8 -> 6 -> 4 -> 0 -> 2
                            if (isVisible(x, y - 1, z)) //Bottom 2
                                _addFace(8);  //FINE
                            if (isVisible(x, y, z - 1)) //Right 4
                                _addFace(0);//_addFace(0);  //FINE					
                            if (isVisible(x, y, z + 1)) //Left 6
                                _addFace(4);//_addFace(4);  //FINE				
                            if (isVisible(x - 1, y, z)) //Back 8
                                _addFace(6);//_addFace(6); //FINE						
                            if (isVisible(x + 1, y, z)) //Front 10
                                _addFace(2); //FINE
                            if (isVisible(x, y + 1, z))
                                _addFace(10);
                }
            }
        }

        glBufferData(GL_ARRAY_BUFFER, chunk->info_count * sizeof(*draw_buffer), draw_buffer, GL_STATIC_DRAW);

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 28, 0);
        glEnableVertexAttribArray(0);

        glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 28, 12);
        glEnableVertexAttribArray(1);

        glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 28, 20);
        glEnableVertexAttribArray(2);

        chunk->needsMeshed = false;
        printf("Uploaded %i verticies\n", chunk->info_count);
    }
}

void game_t::render() {
    drawen_buffers = 0;
    currentWorld->chunks->render();
    printf("Drew %i buffers\n", drawen_buffers);
    /*
    std::vector<chunk_t*> *chunks = &currentWorld->chunks;
    //fprintf(stderr, "Rendering");
    for (auto chunk : *chunks) {
        if (!chunk->loaded) {
            fprintf(stderr, "Not yet loaded\n");
            continue;
        }
        if (chunk->needsMeshed) {
            fprintf(stderr, "Meshing\n");
            mesh(chunk);
        }

        //glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glBindVertexArray(chunk->vao);
        //glBindTexture(GL_TEXTURE_2D, textureId);

        glm::mat4 model = glm::translate(glm::mat4(1.0), glm::vec3(chunk->position.x,chunk->position.y,chunk->position.z));
        glUniformMatrix4fv(uniModel, 1, GL_FALSE, glm::value_ptr(model));

        // Draw the quad
        //glBegin(GL_LINES);
        //printf("Drawing %i verticies\n", chunk->info_count);
        glDrawArrays(GL_TRIANGLES, 0, chunk->info_count);
        //glEnd();
        // Unbind the VAO
        glBindVertexArray(0);

    } 
    */ 
}

inline std::string* readShaderFile(const char* filePath) {
    std::ifstream shaderFile(filePath);
    std::stringstream shaderStream;
    shaderStream << shaderFile.rdbuf();
    std::string *ret = new std::string(shaderStream.str());
    shaderFile.close();
    return ret;
}

bool compileShader(GLuint shaderId, const char* shaderSource) {
    glShaderSource(shaderId, 1, &shaderSource, nullptr);
    glCompileShader(shaderId);

    GLint compileStatus;
    glGetShaderiv(shaderId, GL_COMPILE_STATUS, &compileStatus);
    if (compileStatus == GL_FALSE) {
        GLint logLength;
        glGetShaderiv(shaderId, GL_INFO_LOG_LENGTH, &logLength);
        std::string errorMessage(logLength, ' ');
        glGetShaderInfoLog(shaderId, logLength, nullptr, &errorMessage[0]);
        std::cout << "Shader compilation error: " << errorMessage << std::endl;
        return false;
    }

    return true;
}

bool linkProgram(GLuint programId, GLuint vertexShaderId, GLuint fragmentShaderId) {
    glAttachShader(programId, vertexShaderId);
    glAttachShader(programId, fragmentShaderId);
    glLinkProgram(programId);

    GLint linkStatus;
    glGetProgramiv(programId, GL_LINK_STATUS, &linkStatus);
    if (linkStatus == GL_FALSE) {
        GLint logLength;
        glGetProgramiv(programId, GL_INFO_LOG_LENGTH, &logLength);
        std::string errorMessage(logLength, ' ');
        glGetProgramInfoLog(programId, logLength, nullptr, &errorMessage[0]);
        std::cout << "Program linking error: " << errorMessage << std::endl;
        return false;
    }

    glDetachShader(programId, vertexShaderId);
    glDetachShader(programId, fragmentShaderId);

    return true;
}

int figureOutShaders() {
        // Load and compile vertex shader
    std::string *vertexShaderSource =(readShaderFile("vertex_shader.glsl"));
    const char* vertexShaderSourcePtr = vertexShaderSource->c_str();
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    if (!compileShader(vertexShader, vertexShaderSourcePtr)) {
        glfwTerminate();
        return -1;
    }

    // Load and compile fragment shader
    std::string *fragmentShaderSource = (readShaderFile("fragment_shader.glsl"));
    const char* fragmentShaderSourcePtr = fragmentShaderSource->c_str();
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    if (!compileShader(fragmentShader, fragmentShaderSourcePtr)) {
        glfwTerminate();
        return -1;
    }

    // Create shader program and link shaders
    shaderProgram = glCreateProgram();
    if (!linkProgram(shaderProgram, vertexShader, fragmentShader)) {
        glfwTerminate();
        return -1;
    }

    // Use the shader program
    glUseProgram(shaderProgram);

	uniModel = glGetUniformLocation(shaderProgram, "model");
	
	
	uniView = glGetUniformLocation(shaderProgram, "view");
	
	
	/*glm::mat4*/ proj = glm::perspective(glm::radians(90.0f), 1440.0f / 900.0f, 0.1f, 500.0f);
	//frust.SetFrustum(90.0f, 1440.0f / 900.0f, 0.1f, 500.0f);
    GLint uniProj = glGetUniformLocation(shaderProgram, "proj");
    glUniformMatrix4fv(uniProj, 1, GL_FALSE, glm::value_ptr(proj));
	
    //delete vertexShaderSource;
    //delete vertexShaderSourcePtr;
    //delete fragmentShaderSource;
    //delete fragmentShaderSourcePtr;
}

int main() {
    // Initialize GLFW
    if (!glfwInit()) {
        std::cout << "Failed to initialize GLFW" << std::endl;
        return -1;
    }

		view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);
    // Set error callback
    glfwSetErrorCallback(errorCallback);

    // Create a window
    GLFWwindow* window = glfwCreateWindow(800, 600, "OpenGL Window", NULL, NULL);
    if (!window) {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }

    // Make the window's context current
    glfwMakeContextCurrent(window);

    // Initialize GLEW
    if (glewInit() != GLEW_OK) {
        std::cout << "Failed to initialize GLEW" << std::endl;
        glfwTerminate();
        return -1;
    }

    // Load texture using stb_image
    int width, height, channels;
    unsigned char* image = stbi_load("terrain.png", &width, &height, &channels, 0);
    if (!image) {
        std::cout << "Failed to load texture" << std::endl;
        glfwTerminate();
        return -1;
    }

    // Print texture information
    std::cout << "Loaded texture: width = " << width << ", height = " << height << ", channels = " << channels << std::endl;

    if (figureOutShaders()) {
        std::cout << "Can't into shaders" << std::endl;
        return -1;
    }

        // Load texture into OpenGL
    GLuint textureId;
    glGenTextures(1, &textureId);
    glBindTexture(GL_TEXTURE_2D, textureId);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, image);
    glGenerateMipmap(GL_TEXTURE_2D);

    // Set texture wrapping and filtering options
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);


    // Set the texture unit for the shader
    GLint textureSamplerLocation = glGetUniformLocation(shaderProgram, "textureSampler");
    glUniform1i(textureSamplerLocation, 0);


    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);	
	glfwSetKeyCallback(window, key_frame);
	glfwSetCursorPosCallback(window, mouse_callback);
	//glfwSetMouseButtonCallback(window, mouse_button_callback);

	glfwSwapInterval(1);

    start_game();

	static const GLfloat g_vertex_buffer_data[] = {
    -1.0f,-1.0f,-1.0f, // triangle 1 : begin
    -1.0f,-1.0f, 1.0f,
    -1.0f, 1.0f, 1.0f, // triangle 1 : end
    1.0f, 1.0f,-1.0f, // triangle 2 : begin
    -1.0f,-1.0f,-1.0f,
    -1.0f, 1.0f,-1.0f, // triangle 2 : end
    1.0f,-1.0f, 1.0f,
    -1.0f,-1.0f,-1.0f,
    1.0f,-1.0f,-1.0f,
    1.0f, 1.0f,-1.0f,
    1.0f,-1.0f,-1.0f,
    -1.0f,-1.0f,-1.0f,
    -1.0f,-1.0f,-1.0f,
    -1.0f, 1.0f, 1.0f,
    -1.0f, 1.0f,-1.0f,
    1.0f,-1.0f, 1.0f,
    -1.0f,-1.0f, 1.0f,
    -1.0f,-1.0f,-1.0f,
    -1.0f, 1.0f, 1.0f,
    -1.0f,-1.0f, 1.0f,
    1.0f,-1.0f, 1.0f,
    1.0f, 1.0f, 1.0f,
    1.0f,-1.0f,-1.0f,
    1.0f, 1.0f,-1.0f,
    1.0f,-1.0f,-1.0f,
    1.0f, 1.0f, 1.0f,
    1.0f,-1.0f, 1.0f,
    1.0f, 1.0f, 1.0f,
    1.0f, 1.0f,-1.0f,
    -1.0f, 1.0f,-1.0f,
    1.0f, 1.0f, 1.0f,
    -1.0f, 1.0f,-1.0f,
    -1.0f, 1.0f, 1.0f,
    1.0f, 1.0f, 1.0f,
    -1.0f, 1.0f, 1.0f,
    1.0f,-1.0f, 1.0f
};
    
    GLuint vbo_test;
    glGenBuffers(1, &vbo_test);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_test);
    glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data), g_vertex_buffer_data, GL_STATIC_DRAW);

    // Main loop
    while (!glfwWindowShouldClose(window)) {
        // Render here
		float currentFrame = glfwGetTime();
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;  
        key_callback(window);
		view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);
		glUniformMatrix4fv(uniView, 1, GL_FALSE, glm::value_ptr(view));
        // Swap front and back buffers

        glm::mat4 model = glm::translate(glm::mat4(1.0), glm::vec3(0,0,0));
        glUniformMatrix4fv(uniModel, 1, GL_FALSE, glm::value_ptr(model));
        const glm::vec3 sunLight(0.5f, 0.8f, 1.0f);
        float slSin = sin(0.5) * 0.45f + 0.6f;
        glClearColor(sunLight.x * slSin, sunLight.y * slSin, sunLight.z * slSin, 0.5f);
        
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glEnable(GL_DEPTH_TEST);
		glEnable(GL_POLYGON_OFFSET_FILL);
		
		glEnable(GL_ALPHA_TEST);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		
        
         glUseProgram(shaderProgram);

        printf("\rX: %f, Y: %f, Z: %f", cameraPos.x, cameraPos.y, cameraPos.z);
        
        glBindBuffer(GL_ARRAY_BUFFER, vbo_test);
        glDrawArrays(GL_TRIANGLES, 0, 12*3);

        // Render a simple triangle
        glBegin(GL_TRIANGLES);
            glColor3f(1.0f, 0.0f, 0.0f); // Red
            glVertex3f(-0.6f, -0.4f, 0.0f);
            glColor3f(0.0f, 1.0f, 0.0f); // Green
            glVertex3f(0.6f, -0.4f, 0.0f);
            glColor3f(0.0f, 0.0f, 1.0f); // Blue
            glVertex3f(0.0f, 0.6f, 0.0f);
        glEnd();
        


        currentGame.render();


        // Poll for and process events
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // Free resources
    stbi_image_free(image);

    // Terminate GLFW
    glfwTerminate();

    return 0;
}

void key_frame(GLFWwindow *window) {	
	int cond = GLFW_PRESS;
	if (glfwGetKey(window, GLFW_KEY_1) == cond) {
		selectedBlock--;
		if (selectedBlock < 0)
			selectedBlock = 1;
	}
	if (glfwGetKey(window, GLFW_KEY_2) == cond) {
		selectedBlock = ++selectedBlock % 1;
	}	
	if(glfwGetKey(window, GLFW_KEY_F) == cond)
		glfwSetWindowMonitor(window, monitor, 0, 0, 1440, 900, 0);
	if (glfwGetKey(window, GLFW_KEY_KP_DECIMAL) == cond)
		showTriangles = !showTriangles;
	if (glfwGetKey(window, GLFW_KEY_KP_0) == cond) {
		doDebug = !doDebug;
		/*
        for (int x = 0; x < 3; x++) {
			for (int y = 0; y < 3; y++) {
				worldObj->columns[x][y].vertsChanged = true;
			
			}
		}
        */
	}
}	

void key_callback(GLFWwindow* window) {
	int cond = GLFW_PRESS;
    //puts("Key pressed");
	
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == cond)
        glfwSetWindowShouldClose(window, true);
	float mult = 2.0f;
	if (glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == cond)
		mult = 10.0f;
	
	glm::vec3 cameraFront2(cameraFront.x,0.0f,cameraFront.z);
	
	float cameraSpeed = 5.0f * deltaTime;
	if (glfwGetKey(window, GLFW_KEY_W) == cond)
        cameraPos += cameraSpeed * glm::normalize(cameraFront2) * mult;
    if (glfwGetKey(window, GLFW_KEY_S) == cond)
        cameraPos -= cameraSpeed * glm::normalize(cameraFront2) * mult;
    if (glfwGetKey(window, GLFW_KEY_A) == cond)
        cameraPos -= glm::normalize(glm::cross(cameraFront2, cameraUp)) * cameraSpeed * mult;
    if (glfwGetKey(window, GLFW_KEY_D) == cond)
        cameraPos += glm::normalize(glm::cross(cameraFront2, cameraUp)) * cameraSpeed * mult;
	if (glfwGetKey(window, GLFW_KEY_SPACE) == cond)
		cameraPos.y += cameraSpeed * mult;
	if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == cond)
		cameraPos.y -= cameraSpeed * mult;
	
	//worldObj->updateColumnsAroundPlayer({cameraPos.x, cameraPos.y, cameraPos.z});
	
	if (glfwGetKey(window, GLFW_KEY_Q) == cond) {
		//worldObj->updateColumnsAroundPlayer(worldObj->player.vCamera, true);
	}		
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
	//blockComplete block;
	//printf("block: %f, %f, %f; prev: %f, %f, %f\r\n", worldObj->player.blockLookingAt.x, worldObj->player.blockLookingAt.y, worldObj->player.blockLookingAt.z, worldObj->player.prevBlockLookingAt.x, worldObj->player.prevBlockLookingAt.y, worldObj->player.prevBlockLookingAt.z);
}

void mouse_callback(GLFWwindow* window, double xpos, double ypos) {
    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top
    lastX = xpos;
    lastY = ypos;

    float sensitivity = 0.10f; // change this value to your liking
    xoffset *= sensitivity;
    yoffset *= sensitivity;

    yaw += xoffset;
    pitch += yoffset;

    // make sure that when pitch is out of bounds, screen doesn't get flipped
    if (pitch > 89.9f)
        pitch = 89.9f;
    if (pitch < -89.9f)
        pitch = -89.9f;

    glm::vec3 front;
    front.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
    front.y = sin(glm::radians(pitch));
    front.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
    cameraFront = glm::normalize(front);
}