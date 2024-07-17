#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "include/stb_image.h"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/ext.hpp>
#include <glm/gtx/string_cast.hpp>

#include "include/Shader.h"
#include "include/Camera.h"
#include "include/Mesh.h"
//#include "include/Model.h"

#include <iostream>
#include <string>
#include <cmath>

#include "Simulation.h"
#include "particle.h"

#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>

// settings
const unsigned int SCR_WIDTH = 1000;
const unsigned int SCR_HEIGHT = 1000;

// camera
Camera camera(glm::vec3(5.0f, 5.0f, 15.0f));
float lastX = (float)SCR_WIDTH / 2.0;
float lastY = (float)SCR_HEIGHT / 2.0;
bool firstMouse = true;

// timing
float deltaTime = 0.0f;
float lastFrame = 0.0f;

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void processInput(GLFWwindow* window);

const int FPS = 30;
std::string str_cmd = "ffmpeg -r " + std::to_string(FPS) + " -f rawvideo -pix_fmt rgba -s "
+ std::to_string(SCR_WIDTH) + "x" + std::to_string(SCR_HEIGHT)
+ " -i - -threads 0 -preset fast -y -pix_fmt yuv420p -crf 21 -vf vflip ./sand.mp4";
const char* cmd = str_cmd.c_str();
FILE* ffmpeg = _popen(cmd, "wb");
int* buffer = new int[SCR_WIDTH * SCR_HEIGHT];

std::shared_ptr<Mesh> get_sphere(float r, size_t slices, size_t stacks) {
    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;

    const float kPi = std::atan(1.0f) * 4;

    float phi_step = kPi * 2 / slices;
    float theta_step = kPi / stacks;

    for (size_t vi = 0; vi < stacks; vi++) {  // vertical loop
        float theta = vi * theta_step;
        float z = r * cosf(theta);
        for (size_t hi = 0; hi < slices; hi++) {  // horizontal loop
            float phi = hi * phi_step;
            glm::vec3 p(r * cosf(phi) * sinf(theta), r * sinf(phi) * sinf(theta), z);
            Vertex vertex;
            vertex.Position = p;
            vertex.Normal = glm::normalize(p);
            vertices.push_back(vertex);
        }
    }

    for (size_t vi = 0; vi < stacks; vi++)
        for (size_t hi = 0; hi < slices; hi++) {
            auto t1 = (unsigned int)(vi * slices + hi);
            auto t2 = (unsigned int)(vi * slices + hi + 1);
            auto t3 = (unsigned int)((vi + 1) * slices + hi + 1);
            auto t4 = (unsigned int)((vi + 1) * slices + hi);
            indices.insert(indices.end(), { t1, t2, t3 });
            indices.insert(indices.end(), { t1, t3, t4 });
        }

    std::shared_ptr<Mesh> sphere_mesh = std::make_unique<Mesh>(vertices, indices);

    return sphere_mesh;
}

int main()
{
    // glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // glfw window creation
    // --------------------
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "Water Simulation 2D", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);

    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    // configure global opengl state
    // -----------------------------
    glEnable(GL_DEPTH_TEST);

    // build and compile shaders
    // -------------------------
    Shader shader("./shader/shader.vs", "./shader/shader.fs");
    Shader bordershader("./shader/shader_ez.vs", "./shader/shader_ez.fs");

    // generate a list of 100 quad locations/translation-vectors
    // ---------------------------------------------------------

    //auto sphere = get_square();
    auto sphere = get_sphere(0.1f, 20, 20);

    //==========================================================================================
    //glm::vec2 x0(1.5f, 1.5f, 1.5f);
    glm::vec2 v0(0.0f, 0.0f);
    glm::vec2 v1(-30.0f, -20.0f);
    glm::mat2 B0(0.0f);
    //float mass = 0.0005f;
    //float V0 = 1.0f;

    float dt = 0.001f;

    //std::vector<Water> outParticles;
    std::vector<Sand> outParticles;

    //for (int p = 0; p < 8; p++)
    //{
    //    float r = ((double)rand() / (RAND_MAX));			// random number

    //    int row = p / 5;
    //    int col = p % 5;

    //    //glm::vec2 pos (100.0f, (double)(96.0f - 0.5f * (float)p - r));		// new positions
    //    //glm::vec2 pos(100.0f, 96.0f - 0.5*p - r);
    //    glm::vec2 pos(95.0f, 2.0f + 0.5f * (float)p);
    //    Water the_water(pos, v0, B0, mass, V0);
    //    outParticles.push_back(the_water);
    //}

    //Sand Initialization
    srand((unsigned)time(NULL));
    for (int p = 0; p < 10; p++)
    {
        for (int i = 0; i < 40; i++) {
            float r1 = ((float)rand() / (RAND_MAX));			// random number
            float r2 = ((float)rand() / (RAND_MAX));			// random number
            glm::vec2 pos(40.0f + 0.5f * (float)i + 0.5*r1, 50.0f + 2.5f * (float)p + 2.5 * r2);
            float mass = 32.0f;
            float V0 = 2.0f;
            Sand the_water(pos, v0, B0, mass, V0);
            outParticles.push_back(the_water);
        }
    }

    Simulation simulation(outParticles,dt);

    glm::mat4 projection = glm::perspective(glm::radians(45.0f), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
    glm::mat4 view = camera.GetViewMatrix();

    //==========================================================================================

    Eigen::Matrix2f u, v, cov, r;
    Eigen::Vector2f w;
    Eigen::JacobiSVD<Eigen::Matrix2f> svd;
    glm::mat2 A(1.0f, 0.5f, 0.5f, 1.0f);
    //cov << 1.0f, 0.5f, 0.5f, 1.0f;
    cov << A[0][0], A[0][1], A[1][0], A[1][1];
    svd.compute(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
    u = svd.matrixU(); 
    v = svd.matrixV();
    w = svd.singularValues();
    Eigen::Matrix2f Eps;
    //Eps << w(0), 0.0f, 0.0f, w(1);
    //Eigen::Matrix2f e = Eps.inverse();
    //Eigen::Matrix2f f = Eps.log();
    //float g= f.sum();
    //std::cout << u << endl;
    //std::cout << v << endl;
    //std::cout << w << endl;
    //int a = 1;

    // render loop
    // -----------
    int count = 0;
    while (!glfwWindowShouldClose(window))
    {
        // per-frame time logic
// --------------------
        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        // input
        // -----
        processInput(window);

        // render
        // ------
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        //if (count % 30 == 0) {
        //    for (int i = 0; i < 8; i++) {
        //        float r = ((double)rand() / (RAND_MAX));
        //        glm::vec2 pos(20.0f + (float)i * 0.4f + 0.4f * r, 18.0f);
        //        Water the_water(pos, v0, B0, mass, V0);
        //        simulation.add_particle(the_water);
        //    }
        //
        //srand((unsigned)time(NULL));
        //if (count < 3000 && count % 50 ==0) {
        //    for (int i = 0; i < 30; i++) {
        //        float r = ((double)rand() / (RAND_MAX));
        //        float r2 = ((double)rand() / (RAND_MAX));
        //        float ir = i % 10;
        //        //glm::vec2 pos(95.0f + 0.5 * ir + r, 96.0f - 0.5 * ir - r);
        //        glm::vec2 pos(30.0f + 1.0 * ir + r, 96.0f - r2);
        //        Water the_water(pos, v0, B0, mass, V0);
        //        simulation.add_particle(the_water);
        //    }
        //}

        simulation.forward_one_step();
        if(count % 30 == 0) {
            simulation.render_particles(shader, sphere, projection, view);
            simulation.render_borders(bordershader, projection, view);
            glfwSwapBuffers(window);

            glReadPixels(0, 0, SCR_WIDTH, SCR_HEIGHT, GL_RGBA, GL_UNSIGNED_BYTE, buffer);
            fwrite(buffer, sizeof(int) * SCR_WIDTH * SCR_HEIGHT, 1, ffmpeg);

            glfwPollEvents();
        }
        simulation.reset_grid();

        count++;
        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------

    }
    _pclose(ffmpeg);
    // optional: de-allocate all resources once they've outlived their purpose:
    // ------------------------------------------------------------------------
    //glDeleteVertexArrays(1, &quadVAO);
    //glDeleteBuffers(1, &quadVBO);

    glfwTerminate();
    return 0;
}

// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void processInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        camera.ProcessKeyboard(FORWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        camera.ProcessKeyboard(BACKWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        camera.ProcessKeyboard(LEFT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        camera.ProcessKeyboard(RIGHT, deltaTime);
}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and 
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}

// glfw: whenever the mouse moves, this callback is called
// -------------------------------------------------------
void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
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

    camera.ProcessMouseMovement(xoffset, yoffset);
}

// glfw: whenever the mouse scroll wheel scrolls, this callback is called
// ----------------------------------------------------------------------
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.ProcessMouseScroll(yoffset);
}