#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <fstream>
#include <chrono>
#include <thread>
#include <mutex>
#include <sstream>
#include <unistd.h>
#include <queue>

#include "noise.h"
#include "vertex.h"
#include "glhelper.h"
#include "block.h"

struct dimension_t;
struct game_t;
struct blockaccessor_t;
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
GLuint textureId;
GLuint textureWireframeId;
GLuint fontTextureId;
GLuint textShaderProgram;
bool printDebugInfo = false;
GLFWmonitor* monitor = glfwGetPrimaryMonitor();
void key_callback(GLFWwindow* window);
void key_frame(GLFWwindow* window);
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);

struct game_t {
    struct dimension_t *currentDimension;
    std::vector<vec3d> vectorData;
    void render();
    void mesh(chunk_t *chunk);
    fullblock_t getPlayerLookingAt(vec3d pos, vec3d front);
};

fullblock_t prevBlockLookingAt;

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

static struct _text_draw {
    _text_draw() {
        
    }
    _text_draw(float x, float y, float u, float v) {
        this->x = x;
        this->y = y;
        this->u = u;
        this->v = v;
    }
    float x, y, u, v;
} *text_buffer = 0;

game_t currentGame;
int drawen_buffers;
int triangles;
int meshPerFrame;

const int minimum_octree_size = 2;
const int minimum_block_width = 16;
const int minimum_render_octree = 2;
const int maximum_octree_size = 5;

struct octree_chunk_t : public blockaccessor_t {
    fullblock_t getBlockAtAbsolute(vec3d pos) override {
        //printf(".");


        int width = getBlockWidth(size);
        int halfwidth = width;

        vec3d pos_rel = pos - position;

        fullblock_t nothing_to_return;
        static blockstate_t *nothing_blockstate = new blockstate_t();
        nothing_blockstate->blockId = 1;
        nothing_to_return.block = block_t::stone;
        nothing_to_return.state = nothing_blockstate;

        if (printDebugInfo) {
            fprintf(stderr, "Size: %i, width: %i, pX: %i, pY: %i, pZ: %i, bX: %i, bY: %i, bZ: %i\n", size, width, 
            (int)position.x, (int)position.y, (int)position.z, (int)pos.x, (int)pos.y, (int)pos.z);
        }

        if (!isBlockLoaded(pos)) {
            nothing_to_return.position.x = -15;
            return nothing_to_return;
        }

        if (size > minimum_octree_size) {
            int cwidth = getBlockWidth(size-1);
            int childx = float(pos_rel.x)/float(cwidth);
            int childy = float(pos_rel.y)/float(cwidth);
            int childz = float(pos_rel.z)/float(cwidth);
            if (printDebugInfo) {
                fprintf(stderr, "rX: %f, rY: %f, rZ: %f, cX: %f, cY: %f, cZ: %f\n", pos_rel.x, pos_rel.y, pos_rel.z, childx, childy, childz);
            }
            return children[int((childz * 4) + (childy * 2) + childx)].getBlockAtAbsolute(pos);
            //return nothing_to_return;
        } else {
            if (!point) {
                //printf("{%i,%i,%i} not allocated\n", int(pos.x), int(pos.y), int(pos.z));
                //puts("Nothing here");
                nothing_to_return.position.x = -14;
                return nothing_to_return;
            }

            int x = int(pos_rel.x);
            int y = int(pos_rel.y);
            int z = int(pos_rel.z);
            blockstate_t *state = &point[(y * width * width) + (z * width) + x];
            fullblock_t ret;
            ret.state = state;
            if (state->blockId < block_t::blocks.size())
                ret.block = block_t::blocks.at(state->blockId);
            else
                ret.block = block_t::air;
            ret.position = pos;

            //printf("{%i,%i,%i} not in bounds\n", int(pos.x), int(pos.y), int(pos.z));
            //printf("\n");
            return ret;
        }

        //printf("{%i,%i,%i} we shouldn't be here\n", int(pos.x), int(pos.y), int(pos.z));
        nothing_to_return.position.x = -13;
        return nothing_to_return;
    }

    fullblock_t getBlockAtRelative(vec3d pos) override {
        return getBlockAtAbsolute(pos + position);
    }

    bool isAvailable() override {
        return true;
    }

    bool isBlockLoaded(vec3d pos) override {
        int width = getBlockWidth(size);
        if (pos.x < position.x ||
            pos.x >= position.x + width ||
            pos.y < position.y||
            pos.y >= position.y + width ||
            pos.z < position.z ||
            pos.z >= position.z + width) {
                return false;
            }
        return true;
    }

    octree_chunk_t() {
        has_children = false;
        is_meshable = false;
        needsMeshed = true;
        waitingOnThread = false;
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
        waitingOnThread = false;
        children = 0;
        point = 0;
        vbo = 0;
        vao = 0;
        info_count = 0;
        this->position = _position;

        if (size >= minimum_octree_size) { //If 4x4x4 and above, we hold blocks and can render individually
            glGenBuffers(1, &vbo);
            glGenVertexArrays(1, &vao);
            is_meshable = true;
        }
        if (size == minimum_octree_size) { //If 4x4x4, we hold blocks
            int width = getBlockWidth(size);
            int volume = pow(width, 3);
            point = new blockstate_t[volume];
            memset(point, 0, sizeof(*point) * volume);
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
        if (size > minimum_octree_size) { //All sizes above 4x4x4 create children
            generate_children();
        }
    }

    void generate() {
        int width = getBlockWidth(size);
        for (int y = 0; y < width; y++)
            for (int x = 0; x < width; x++)
                for (int z = 0; z < width; z++) {
                    //if (position.y + y > 0 && position.y + y < 5)
                    //point[(y * width * width) + (z * width) + x] = block_t::stone->getDefaultState();
                }
                    //*getBlockAtRelative(vec3d{x,y,z}).state = block_t::stone->getDefaultState();
//return;
        for (int x = 0; x < width; x++) {
            for (int z = 0; z < width; z++) {
                int ofx = x + position.x;
                int ofz = z + position.z;
                perlin::octaves = 2;
                perlin::persistence = 0.55f;
                //float scale = 0.005f;
                float scale = 0.01f;
                float height = 600.0f;
                float value = perlin::getPerlin(ofx * scale + 16000.0f, ofz * scale + 16000.0f) * height + 10;

                for (int y = 0; y < width; y++) {
                    int ofy = y + position.y;
                    auto sampleSphere = [](double radius, double x, double y, double z, double centerX, double centerY, double centerZ) {		
                        double fx = x - centerX, fy = y - centerY, fz = z - centerZ;
                        return radius * radius > ((fx * fx) + (fy * fy) + (fz * fz));
                    };
                    
                    //if	(sampleSphere(getBlockWidth(maximum_octree_size-1)-1, ofx, ofy, ofz, getBlockWidth(maximum_octree_size-1), getBlockWidth(maximum_octree_size-1), getBlockWidth(maximum_octree_size-1)))
                    //    point[(y * width * width) + (z * width) + x] = block_t::dirt->getDefaultState();
                    //continue;
                    if (value - 2.0f > position.y + y)  
                        point[(y * width * width) + (z * width) + x] = block_t::stone->getDefaultState();
                    else
                        if (value - 1.0f > position.y + y)
                           point[(y * width * width) + (z * width) + x] = block_t::grass->getDefaultState();
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
    bool waitingOnThread;
    int info_count;
    vec3d position;
    static const int info_max = 800000;

    static int getBlockWidth(int size) {
        //float fact = float(minimum_block_width) / pow(2, 2);
        //return int(pow(2,size)) * fact;        
        //return (minimum_block_width * size) / 2.0f;
        return 0.25f*minimum_block_width*(1 << size);
    }

    bool needsToRender() {
        if (!is_meshable || !children || size < minimum_render_octree)
            return false;
        //return true;
        float threshold = getBlockWidth(size) * 4;
        float dist = 
        sqrtf(((cameraPos.x - position.x) * (cameraPos.x - position.x)) +
              ((cameraPos.y - position.y) * (cameraPos.y - position.y)) +
              ((cameraPos.z - position.z) * (cameraPos.z - position.z)));

        //return (dist < pow(2, size-1));
        return dist < threshold;
    }

    void free_children() {
        if (has_children && children) {
            for (int i = 0; i < 8; i++) {
                children[i].free_children();
                children[i].~octree_chunk_t();
            }
        }

        has_children = false;
    }

    void generate_children() {
        //printf("Generating size %i\n", size);
        {
            double childSize = getBlockWidth(size-1);//std::cbrt(pow(8, size-1));
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
        free_children();
        free(children);
        glDeleteBuffers(1, &vbo);
        glDeleteVertexArrays(1, &vao);
    }

    bool independent = false;

    void is_individually_rendering() {
        if (!independent) {
            //puts("Now in range");
            independent = true;

            //for (int i = 0; i < 8; i++)
            //    children[i].is_dependent_meshing();
        }
    }

    void is_dependent_meshing() {
        if (independent) {
            //puts("Now out of range");
            independent = false;
        }
    }

    void mesh(_draw_type *buffer, int *index);

    void mesh_entrypoint() {
        //waitingOnThread = true;

        {
			//std::lock_guard<std::mutex> lk(drawBufferMux);
            if (!draw_buffer) {
                printf("Allocating %i entries\n", info_max);
                draw_buffer = (_draw_type*)malloc(sizeof(_draw_type) * info_max);
            }
            needsMeshed = false;
            info_count = 0;

            //fprintf(stderr, "Meshing\n");            
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
            //waitingOnThread = false;
            printf("Uploaded %i verticies\n", info_count);
        }
    }

    bool out_of_range() {
        int dist = getBlockWidth(size) * 2;
        return abs(cameraPos.x - position.x) > dist ||
               abs(cameraPos.y - position.y) > dist ||
               abs(cameraPos.z - position.z) > dist;
        
        //float dist = 
        //sqrtf(((cameraPos.x - position.x) * (cameraPos.x - position.x)) +
        //      ((cameraPos.y - position.y) * (cameraPos.y - position.y)) +
        //      ((cameraPos.z - position.z) * (cameraPos.z - position.z)));
        //return dist > getBlockWidth(size) * 2.5f;//2.8284f;
    }

    void render() {
        //if (!is_meshable)
        //    return;

        //if (!chunk->loaded) {
        //    fprintf(stderr, "Not yet loaded\n");
        //    continue;
        //}

        { //Test if any children are in the range of their LOD, if they are render to them, ignore ourselves
            bool childRendered = false;
            bool childNeedLOD = false;
            for (int i = 0; i < 8; i++) {
                if (children[i].needsToRender()) {
                    //children[i].render();
                    childRendered = true;
                } else {
                    childNeedLOD = true;
                }
            }

            if (!childNeedLOD && childRendered) {
                for (int i = 0; i < 8; i++)
                    children[i].render();
                is_dependent_meshing();
                return;
            }
        }

        is_individually_rendering(); //Now in range

        if (!needsToRender())
            return;
        
        if (waitingOnThread)
            return;

        if (needsMeshed) {
            if (meshPerFrame < 1)
                mesh_entrypoint();
            if (info_count > 0)
                meshPerFrame++;
            //waitingOnThread = true;
            //workerTask_t *task = new workerTask_t;
			//task->arg1 = this;
			//task->workerTaskFunction = (void(*)(void*))(void(*)(octree_chunk_t*))&octree_chunk_t::mesh_entrypoint;
			//workerTasks.push(task);
            //return;//come back when meshed
        }


        if (info_count == 0)
            return;

        int width = getBlockWidth(size);
            // Calculate the half-width of the cube
        float halfWidth = width;

        // Calculate the coordinates of the eight vertices of the cube
        float xMin = position.x;
        float xMax = position.x + halfWidth;
        float yMin = position.y;
        float yMax = position.y + halfWidth;
        float zMin = position.z;
        float zMax = position.z + halfWidth;

        // Define the edges of the cube
        // Each line connects two vertices
        std::vector<std::pair<glm::vec3, glm::vec3>> edges = {
            // Bottom face
            {glm::vec3(xMin, yMin, zMin), glm::vec3(xMax, yMin, zMin)},
            {glm::vec3(xMin, yMin, zMin), glm::vec3(xMin, yMin, zMax)},
            {glm::vec3(xMax, yMin, zMin), glm::vec3(xMax, yMin, zMax)},
            {glm::vec3(xMin, yMin, zMax), glm::vec3(xMax, yMin, zMax)},

            // Top face
            {glm::vec3(xMin, yMax, zMin), glm::vec3(xMax, yMax, zMin)},
            {glm::vec3(xMin, yMax, zMin), glm::vec3(xMin, yMax, zMax)},
            {glm::vec3(xMax, yMax, zMin), glm::vec3(xMax, yMax, zMax)},
            {glm::vec3(xMin, yMax, zMax), glm::vec3(xMax, yMax, zMax)},

            // Vertical edges
            {glm::vec3(xMin, yMin, zMin), glm::vec3(xMin, yMax, zMin)},
            {glm::vec3(xMax, yMin, zMin), glm::vec3(xMax, yMax, zMin)},
            {glm::vec3(xMin, yMin, zMax), glm::vec3(xMin, yMax, zMax)},
            {glm::vec3(xMax, yMin, zMax), glm::vec3(xMax, yMax, zMax)},
        };

        glUseProgram(shaderProgram);
        glBindTexture(GL_TEXTURE_2D, textureWireframeId);
        // Use the edges to draw the wireframe cube
        glBegin(GL_LINES);
        for (const auto& edge : edges) {
            const glm::vec3& v1 = edge.first;
            const glm::vec3& v2 = edge.second;
            glVertex3f(v1.x, v1.y, v1.z);
            glVertex3f(v2.x, v2.y, v2.z);
        }
        glEnd();
/*
        glBegin(GL_LINES);
        glBindTexture(GL_TEXTURE_2D, textureWireframeId);
        // Draw the cube edges
        glVertex3f(position.x, position.y, position.z);
        glVertex3f(position.x, position.z, position.z + width);
        
        glVertex3f(position.x, position.y, position.z + width);
        glVertex3f(position.x, position.y + width, position.z + width);
        
        // ... Continue adding lines for each edge of the cube
        glEnd();
        */

        glBindTexture(GL_TEXTURE_2D, textureId);
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
        triangles += info_count;
        //glEnd();
        // Unbind the VAO
        glBindVertexArray(0);

    }
};

struct chunk_cache_t /*: public blockaccessor_t*/ {
    std::vector<octree_chunk_t*> cache;
    std::mutex cacheMutex;

    bool is_cached(vec3d origin) {
        {
            std::unique_lock<std::mutex> lock(cacheMutex);
            for (octree_chunk_t* chunk : cache)
                if (chunk->position.equals(origin, 0.1f))
                    return true;
            return false;
        }
    }

    octree_chunk_t *get_cached(vec3d origin) {
        {
            std::unique_lock<std::mutex> lock(cacheMutex);
            for (octree_chunk_t* chunk : cache)
                if (chunk->position.equals(origin, 0.1f))
                    return chunk;
        }
        return 0;
    }

    void add_cached(octree_chunk_t* chunk) {
        {
            std::unique_lock<std::mutex> lock(cacheMutex);
            cache.push_back(chunk);
        }
    }

    octree_chunk_t *get_chunk(vec3d origin, bool getCachedOnly, bool allowGeneration) {
        //first check vector
        if (is_cached(origin))
            return get_cached(origin);

        if (getCachedOnly)
            return 0;

        if (!allowGeneration)
            return 0;

        printf("Created new chunk x:%.1f y:%.1f z:%.1f\n", origin.x, origin.y, origin.z);
        octree_chunk_t *new_chunk = new octree_chunk_t(maximum_octree_size, origin);
        add_cached(new_chunk);
        return new_chunk;
    }

    void free_chunk(vec3d origin) {
        //write to disk or whatever
    }
};

struct dimension_t : public blockaccessor_t {
    //std::vector<chunk_t*> chunks;
    const int slot_count = 27;
    //octree_chunk_t *chunks = new octree_chunk_t(maximum_octree_size, vec3d{0,0,0});
    //octree_chunk_t *chunks[27];

    chunk_cache_t dimension_chunk_cache;

    void set_chunks() {
        int width = octree_chunk_t::getBlockWidth(maximum_octree_size);
        int wwidth = width * 3;
        int half = wwidth / 2;

        //for (int i = 0; i < slot_count; i++) {
        //    if (chunks[i] && chunks[i]->out_of_range()) {
        //        puts("chunk out of range");
        //        chunks[i]->~octree_chunk_t();
        //        chunks[i] = 0;
        //    }
        //}
                        
        //int playerrx0 = ((int)cameraPos.x - half) / width;
        //int playerry0 = ((int)cameraPos.y - half) / width;
        //int playerrz0 = ((int)cameraPos.z - half) / width;
        int playerrx0 = ((int)cameraPos.x) / width;
        int playerry0 = ((int)cameraPos.y) / width;
        int playerrz0 = ((int)cameraPos.z) / width;

        { //Test around the player and free unused chunks
            for (octree_chunk_t *chunk : dimension_chunk_cache.cache) { //may need to lock later
                bool has = false;
                for (int xr = 0; xr < 3; xr++)
                    for (int yr = 0; yr < 3; yr++)
                        for (int zr = 0; zr < 3; zr++) {
                            int x = ((-((xr) & 1) * 2)+1) * ((xr+1)/2);
                            int y = ((-((yr) & 1) * 2)+1) * ((yr+1)/2);
                            int z = ((-((zr) & 1) * 2)+1) * ((zr+1)/2);
                            int testx = (playerrx0 + x) * width;
                            int testy = (playerry0 + y) * width;
                            int testz = (playerrz0 + z) * width;

                            vec3d test_origin = vec3d(testx, testy, testz);

                            if (dimension_chunk_cache.is_cached(test_origin)) {
                                has = true;
                                goto done;
                                break;
                            }
                        }
                done:;
                if (!has) {
                    printf("Does not have chunk cached. Free this chunk\n");
                }
            }
        }

        { //Test around the player and create chunks if they are not cached
            for (int xr = 0; xr < 3; xr++) {
                for (int yr = 0; yr < 3; yr++) {
                    for (int zr = 0; zr < 3; zr++) {
                        int x = ((-((xr) & 1) * 2)+1) * ((xr+1)/2);
                        int y = ((-((yr) & 1) * 2)+1) * ((yr+1)/2);
                        int z = ((-((zr) & 1) * 2)+1) * ((zr+1)/2);
                        int testx = (playerrx0 + x) * width;
                        int testy = (playerry0 + y) * width;
                        int testz = (playerrz0 + z) * width;

                        vec3d test_origin = vec3d(testx, testy, testz);

                        bool has = dimension_chunk_cache.is_cached(test_origin);

                        if (has)
                            continue;

                        if (dimension_chunk_cache.get_chunk(test_origin, false, true)) {
                            puts("Allocated new chunk");
                        } else {
                            puts("Failed to allocate new chunk, return address 0");
                        }
                    }
                }
            }
        }
    }

    dimension_t() {
        int width = octree_chunk_t::getBlockWidth(maximum_octree_size);
        int wwidth = width * 3;
        int half = wwidth / 2;
        int playerx0 = ((int)cameraPos.x % wwidth) - half;
        int playery0 = ((int)cameraPos.y % wwidth) - half;
        int playerz0 = ((int)cameraPos.z % wwidth) - half;
        /*
        for (int x = 0; x < 3; x++) {
            for (int y = 0; y < 3; y++) {
                for (int z = 0; z < 3; z++) {
                    chunks[(x * 9) + (y * 3) + z] = new octree_chunk_t(maximum_octree_size, vec3d{
                        playerx0 + (x * width),
                        playery0 + (y * width),
                        playerz0 + (z * width)
                    });
                }
            }
        }
        */
    }

    fullblock_t getBlockAtAbsolute(vec3d pos) override {
        return getBlockAtRelative(pos);
    }

    fullblock_t getBlockAtRelative(vec3d pos) override {
        //for (chunk_t* chunk : chunks) {
        //    if (chunk->isBlockLoaded(pos))
        //        return chunk->getBlockAtAbsolute(pos);
        //}
        
        /*
        for (int i = 0; i < slot_count; i++)
            if (chunks[i] && chunks[i]->isBlockLoaded(pos))
                return chunks[i]->getBlockAtAbsolute(pos);
        */
        for (octree_chunk_t *chunk : dimension_chunk_cache.cache)
            if (chunk->isBlockLoaded(pos))
                return chunk->getBlockAtAbsolute(pos);

        //default
        fullblock_t nothing_fb;
        nothing_fb.block = block_t::stone;
        nothing_fb.state = &nothing;
        return nothing_fb;
        //return chunks[0]->getBlockAtAbsolute(pos);

        //return chunks->getBlockAtRelative(pos);
    }

    bool isBlockLoaded(vec3d pos) override {
        /*
        for (int i = 0; i < slot_count; i++)
            if (chunks[i] && chunks[i]->isBlockLoaded(pos))
                return true;
        */
        for (octree_chunk_t *chunk : dimension_chunk_cache.cache)
            if (chunk->isBlockLoaded(pos))
                return true;
        return false;
        //return chunks->isBlockLoaded(pos);
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
                                float textureAtlas[] = {1*factor,0*factor,2*factor,1*factor};
                                //SHOWS ONLY ONES INDIVIDUALLY RENDERED BECAUSE OF APOSITION SHADER
                                draw_buffer[*index].a.x = positionCube3[(l * 9) + (l1 * 3) + 0] + ofx;// + (float(ofx) / 128.0f);
                                draw_buffer[*index].a.y = positionCube3[(l * 9) + (l1 * 3) + 1] + ofy;// + (float(ofy) / 128.0f);
                                draw_buffer[*index].a.z = positionCube3[(l * 9) + (l1 * 3) + 2] + ofz;// + (float(ofz) / 128.0f);
                                
                                draw_buffer[*index].a.u = textureAtlas[map[l % 2][l1][0]];
                                draw_buffer[*index].a.v = textureAtlas[map[l % 2][l1][1]];
                                
                                
                                draw_buffer[*index].a.i = 1.0f;
                                draw_buffer[*index].a.a = 1.0f; //world brightness
                                
                                (*index)++;
                            }
                        }

        //If we have children that may use their vbos, we either add them to our mesh or let them handle it
        if (size > minimum_octree_size) {
            bool childNotMeshed = true;
            if (has_children) {
                for (int i = 0; i < 8; i++) {
                    octree_chunk_t* child = &children[i];
                    //if (!child->needsToRender()) {
                        child->is_dependent_meshing();
                       // child->mesh(buffer, index); //Nothing renders if this is commented out
                    //} else {
                    //    childNotMeshed = false;
                    //}
                }
            }

            //if (!childNotMeshed)
           //     return; //If no child contributed to the mesh, we need to draw ourselves
        }

        //if (size > minimum_octree_size)
        //    return;

        //We have blocks that can be meshed
        //puts("Mesh");
        //printf("%i %i %i\n", ofx, ofy, ofz);

        int width = getBlockWidth(size);
        int skip = size - minimum_octree_size;
        int block_scale = skip * 3;
        if (skip == 1)
            block_scale = 1;
        for (int x = 0; x < width; x+=skip) {
            for (int y = 0; y < width; y+=skip) {
                for (int z = 0; z < width; z+=skip) {
                    
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
                        fullblock_t _fb = currentGame.currentDimension->getBlockAtAbsolute(vec3d{_x,_y,_z});
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
                        
                        for (int l = _ofs * 2; l < _ofs * 2 + 2; l++) {
                            for (int l1 = 0; l1 < 3; l1++) {			
                                if ((*index) + 2 > info_max) {
                                    puts("Could not resize buffers");
                                    //exit(1);b
                                    printf("%i/%i %p\n", *index, info_max, index);
                                }
                            
                                float factor = (1.0f/16.0f);
                                //float textureAtlas[] = {7*factor,0*factor,8*factor,1*factor};
                                float textureAtlas[] = {fb.block->textureAtlas[0]*factor,
                                                        fb.block->textureAtlas[1]*factor,
                                                        fb.block->textureAtlas[2]*factor,
                                                        fb.block->textureAtlas[3]*factor};

                                draw_buffer[*index].a.x = positionCube3[(l * 9) + (l1 * 3) + 0] * block_scale + float(x + ofx);
                                draw_buffer[*index].a.y = positionCube3[(l * 9) + (l1 * 3) + 1] * block_scale + float(y + ofy);
                                draw_buffer[*index].a.z = positionCube3[(l * 9) + (l1 * 3) + 2] * block_scale + float(z + ofz);
                                
                                draw_buffer[*index].a.u = fb.block->getTexture(l,l1,0) * factor;//textureAtlas[map[l % 2][l1][0]];
                                draw_buffer[*index].a.v = fb.block->getTexture(l,l1,1) * factor;//textureAtlas[map[l % 2][l1][1]];
                                
                                
                                draw_buffer[*index].a.i = ((float)fb.block->getFaceBrightness(_ofs,fb))/255.0f; //block brightness
                                draw_buffer[*index].a.a = ((float)fb.block->getFaceBrightness(_ofs,fb))/255.0f; //world brightness
                                
                                (*index)++;
                            }
                        }
                        
                    };
                    
                    //Old order  10 -> 8 -> 0 -> 4 -> 6 -> 2
                    //New order? 10 -> 8 -> 6 -> 4 -> 0 -> 2
                            if (isVisible(x, y - 1, z)) //Bottom 2
                                _addFace(BOTTOM_FACE);  //FINE
                            if (isVisible(x, y, z - 1)) //Right 4
                                _addFace(RIGHT_FACE);//_addFace(0);  //FINE					
                            if (isVisible(x, y, z + 1)) //Left 6
                                _addFace(LEFT_FACE);//_addFace(4);  //FINE				
                            if (isVisible(x - 1, y, z)) //Back 8
                                _addFace(BACK_FACE);//_addFace(6); //FINE						
                            if (isVisible(x + 1, y, z)) //Front 10
                                _addFace(FRONT_FACE); //FINE
                            if (isVisible(x, y + 1, z))
                                _addFace(TOP_FACE);

                if ((*index) + 2 > info_max)
                    goto no_resize;

                }
            }
        }

        no_resize:;

        //size == 2, mesh blocks
}

fullblock_t game_t::getPlayerLookingAt(vec3d origin, vec3d front) {
    fullblock_t block;
    for (float ray = 0.0f; ray < 5.0f; ray += 0.1f) {
        float rx, ry, rz;
        float bx, by, bz;
        //x = sinCFXPI * ray;
        //y = sinCFYPI * ray;
        //z = sinCFZPI * ray;
        rx = front.x * ray;
        ry = front.y * ray;
        rz = front.z * ray;
        
        bx = rx + origin.x - 1;
        by = ry + origin.y - 1;
        bz = rz + origin.z - 1;
        
        //printDebugInfo = true;
        block = currentDimension->getBlockAtAbsolute(vec3d{rx + cameraPos.x, ry + cameraPos.y, rz + cameraPos.z});
        //printDebugInfo = false;
        //printf("\r\t\t\t\t%f %f %f", bx, by, bz);
        
        if (block.state->blockId == 0) {
            prevBlockLookingAt = block;
            if (!(ray + 0.1f < 5.0f)) {
                //worldObj->player.prevBlockLookingAt = { -99999, 99999, 12345 };
            }
            continue;
        }
    }
    return block;
}

struct text_renderer_t {
    std::string string_buffer;
    GLuint vbo, vao; 
    int currentY, currentX;
    float screenX, screenY;
    bool changed;
    int info_count;
    const int vertex_max = 10000;
    

    text_renderer_t() {
        vbo = 0;
        vao = 0;
        changed = true;
        info_count = 0;
        currentX = 0;
        currentY = 0;
        screenX = 0;
        screenY = 0;

        glGenBuffers(1, &vbo);
        glGenVertexArrays(1, &vao);

        if (!text_buffer) {
            text_buffer = new _text_draw[vertex_max];
        }

        //string_buffer = std::string("Hello World!\n2 lines!");
    }

    void set_string(const char *str) {
        string_buffer = std::string(str);
        changed = true;
    }

    void add_character(int character) {
        const float fontWidth = 1.0f / 16.0f;
        const float fontHeight = 1.0f / 16.0f;
        const float screenWidth = 1.0f / 16.0f;
        const float screenHeight = 1.0f / 16.0f;
        float sX = screenX - 1.0f;
        float sY = screenY - (1.0f - screenHeight);
        _text_draw verticies[6] = {
            {sX,sY,0,0},
            {sX+screenWidth,sY,fontWidth,0},
            {sX,sY+screenHeight,0,fontHeight},
            {sX+screenWidth,sY,fontWidth,0},
            {sX+screenWidth,sY+screenHeight,fontWidth,fontHeight},
            {sX,sY+screenHeight,0,fontHeight}
        };
        
        float textX = (character % 16) * fontWidth;
        float textY = (character / 16 % 16) * fontHeight;
        float posX = currentX * screenWidth;
        float posY = currentY * screenHeight;

        for (int i = 0; i < 6; i++) {
            text_buffer[info_count].x = verticies[i].x + posX;
            text_buffer[info_count].y = -verticies[i].y + (screenHeight - posY);
            text_buffer[info_count].u = verticies[i].u + textX;
            text_buffer[info_count].v = verticies[i].v + textY; 

            info_count++;
        }
    }

    void mesh() {
        //fprintf(stderr, "Meshing Text Buffer\n");            
        info_count = 0;
        currentX = 0;
        currentY = 0;
        memset(text_buffer, 0, sizeof(*text_buffer) * vertex_max);

        for (auto ch : string_buffer) {
            if (ch == '\n') {
                currentY++;
                currentX = 0;
                continue;
            }
                
            add_character(ch);
            currentX++;
        }

        changed = false;

        if (info_count == 0)
            return;

        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);

        glBufferData(GL_ARRAY_BUFFER, info_count * sizeof(*text_buffer), text_buffer, GL_STATIC_DRAW);

        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 16, 0);
        glEnableVertexAttribArray(0);

        glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 16, 8);
        glEnableVertexAttribArray(1);

        //printf("Uploaded %i text verticies\n", info_count);
    }

    void render() {
        if (changed)
            mesh();

        glUseProgram(textShaderProgram);
        glBindTexture(GL_TEXTURE_2D, fontTextureId);
        glBindVertexArray(vao);
        glm::mat4 model = glm::translate(glm::mat4(1.0), vec3d{0,0,0});
        glUniformMatrix4fv(uniModel, 1, GL_FALSE, glm::value_ptr(model));
        glDrawArrays(GL_TRIANGLES, 0, info_count);
        glBindVertexArray(0);
    }
} *debug_info;

void start_game() {
    puts("Starting game");
    printf("minimum_block_width:%i\n", minimum_block_width);
    printf("map width: %i\n", octree_chunk_t::getBlockWidth(maximum_octree_size));
    printf("minecraft chunk equivalent: %f\n", octree_chunk_t::getBlockWidth(maximum_octree_size) / 32.0f);
    for (int i = 2; i < 8; i++) {
        printf("%i:%i ", i, octree_chunk_t::getBlockWidth(i));
    }
    printf("\n");
    //exit(1);

    block_t::air = new block_t(0, 0,0,0,0);
    block_t::grass = new block_t(1, 0, 0, 1, 1);
    block_t::stone = new block_t(2, 1, 0, 2, 1);
    block_t::dirt = new block_t(3, 2, 0, 3, 1);

    currentGame.currentDimension = new dimension_t();

    debug_info = new text_renderer_t();
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

void game_t::render() {
    glUseProgram(shaderProgram);
    drawen_buffers = 0;
    meshPerFrame = 0;
    triangles = 0;
    //currentWorld->chunks->render();
    currentDimension->set_chunks();
    //for (int i = 0; i < currentDimension->slot_count; i++)
    //    if (currentDimension->chunks[i])
    //        currentDimension->chunks[i]->render();
    for (octree_chunk_t *chunk : currentDimension->dimension_chunk_cache.cache)
        chunk->render();
    glDisable(GL_CULL_FACE);
    debug_info->render();
    glEnable(GL_CULL_FACE);
    glUseProgram(shaderProgram);
    //printf("Drew %i buffers\n", drawen_buffers);
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

int main() {
    int frameCount = 0;
    double fps;
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;

    start = std::chrono::high_resolution_clock::now();
    // Initialize GLFW
    if (!glfwInit()) {
        std::cout << "Failed to initialize GLFW" << std::endl;
        return -1;
    }

		view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);
    // Set error callback
    glfwSetErrorCallback(errorCallback);

    // Create a window
    GLFWwindow* window = glfwCreateWindow(800, 450, "OpenGL Window", NULL, NULL);
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

    if (!(textShaderProgram = figureOutShaders("text_"))) {
        std::cout << "Can't load text shaders" << std::endl;
        return -1;
    }

    GLuint textuniproj = glGetUniformLocation(textShaderProgram, "projection");
    glm::mat4 text_projection = glm::ortho(0.0f, 800.0f, 0.0f, 600.0f);
    glUniformMatrix4fv(textuniproj, 1, GL_FALSE, glm::value_ptr(text_projection));

    fontTextureId = loadTexture("text.png");

    GLint fontTextureSamplerLocation = glGetUniformLocation(textShaderProgram, "text");
    glUniform1i(fontTextureSamplerLocation, 0);

    if (!(shaderProgram = figureOutShaders(""))) {
        std::cout << "Can't into shaders" << std::endl;
        return -1;
    }

	uniModel = glGetUniformLocation(shaderProgram, "model");
	uniView = glGetUniformLocation(shaderProgram, "view");
	proj = glm::perspective(glm::radians(90.0f), 1920.0f / 1080.0f, 0.1f, 1000.0f);
    GLint uniProj = glGetUniformLocation(shaderProgram, "proj");
    glUniformMatrix4fv(uniProj, 1, GL_FALSE, glm::value_ptr(proj));
	
    textureId = loadTexture("terrain.png");
    textureWireframeId = loadTexture("wireframe_color.png");

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

    glEnable(GL_CULL_FACE);
    glCullFace(GL_FRONT);
    glFrontFace(GL_CW);
	
    // Main loop
    while (!glfwWindowShouldClose(window)) {
        // Render here
		float currentFrame = glfwGetTime();
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;  
        key_callback(window);

        glUseProgram(shaderProgram);
		view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);
		glUniformMatrix4fv(uniView, 1, GL_FALSE, glm::value_ptr(view));
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
		
        currentGame.render();
        

        //printf("\rX: %f, Y: %f, Z: %f", cameraPos.x, cameraPos.y, cameraPos.z);
        frameCount++;
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        if (duration.count() >= 1.0) {
            // Calculate the FPS for the last second
            fps = frameCount / duration.count();

            // Reset the frame count and timer
            frameCount = 0;
            start = std::chrono::high_resolution_clock::now();
        }

        char char_buf[1000];
        fullblock_t lookingAt = currentGame.getPlayerLookingAt(vec3d{cameraPos.x, cameraPos.y, cameraPos.z}, vec3d{cameraFront.x, cameraFront.y, cameraFront.z});
        snprintf(char_buf, 1000, "FPS: %.0f\n<%.2f,%.2f,%.2f>\nvbos: %i chunks: %i\nLooking At: <%i,%i,%i>:%i\nTriangles: %i", fps, cameraPos.x, cameraPos.y, cameraPos.z, drawen_buffers, currentGame.currentDimension->dimension_chunk_cache.cache.size(),
        (int)lookingAt.position.x, (int)lookingAt.position.y, (int)lookingAt.position.z, lookingAt.state->blockId, triangles / 3);
        debug_info->set_string(&char_buf[0]);

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
        




        // Poll for and process events
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // Free resources
    //stbi_image_free(image);

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
    if (glfwGetKey(window, GLFW_KEY_R) == cond) {
        //currentGame.currentWorld->chunks->~octree_chunk_t();
        //currentGame.currentWorld->chunks = new octree_chunk_t(maximum_octree_size, vec3d{
        //cameraPos.x - octree_chunk_t::getBlockWidth(maximum_octree_size-1),
        //cameraPos.y - octree_chunk_t::getBlockWidth(maximum_octree_size-1),
        //cameraPos.z - octree_chunk_t::getBlockWidth(maximum_octree_size-1)});
    }
	if(glfwGetKey(window, GLFW_KEY_F) == cond)
		glfwSetWindowMonitor(window, monitor, 0, 0, 1920, 1080, 0);
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
    fullblock_t lookingAt = currentGame.getPlayerLookingAt(vec3d{cameraPos.x, cameraPos.y, cameraPos.z}, vec3d{cameraFront.x, cameraFront.y, cameraFront.z});
    if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
		if (currentGame.currentDimension->isBlockLoaded(prevBlockLookingAt.position)) {
            *currentGame.currentDimension->getBlockAtAbsolute(prevBlockLookingAt.position).state = block_t::dirt->getDefaultState();
		}
	}
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
		if (currentGame.currentDimension->isBlockLoaded(lookingAt.position)) {
            *currentGame.currentDimension->getBlockAtAbsolute(lookingAt.position).state = block_t::air->getDefaultState();
		}
	}	
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