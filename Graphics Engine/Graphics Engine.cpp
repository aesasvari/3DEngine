#include "olcConsoleGameEngine.h"
#include <strstream>
#include <fstream>
#include <algorithm>
#include <string>

using namespace std;

struct vec3d
{
	float x = 0, y = 0, z = 0, w = 1;
};

struct vec2d
{
	float u = 0, v = 0, w = 1;
};

struct triangle
{
	vec3d p[3];
	vec2d t[3];
	wchar_t sym;
	short col;
};

struct mesh
{
	std::vector<triangle> triangles;

	bool loadObj(string sFilename, bool bHasTexture = false)
	{
		ifstream f(sFilename);
		if (!f.is_open())
			return false;

		// Local cache of verts
		vector<vec3d> verts;
		vector<vec2d> texs;

		while (!f.eof())
		{
			char line[128];
			f.getline(line, 128);

			strstream s;
			s << line;

			char junk;

			if (line[0] == 'v')
			{
				if (line[1] == 't')
				{
					vec2d v;
					s >> junk >> junk >> v.u >> v.v;
					// A little hack for the spyro texture
					//v.u = 1.0f - v.u;
					//v.v = 1.0f - v.v;
					texs.push_back(v);
				}
				else
				{
					vec3d v;
					s >> junk >> v.x >> v.y >> v.z;
					verts.push_back(v);
				}
			}

			if (!bHasTexture)
			{
				if (line[0] == 'f')
				{
					int f[3];
					s >> junk >> f[0] >> f[1] >> f[2];
					triangles.push_back({ verts[f[0] - 1], verts[f[1] - 1], verts[f[2] - 1] });
				}
			}
			else
			{
				if (line[0] == 'f')
				{
					s >> junk;

					string tokens[6];
					int nTokenCount = -1;


					while (!s.eof())
					{
						char c = s.get();
						if (c == ' ' || c == '/')
							nTokenCount++;
						else
							tokens[nTokenCount].append(1, c);
					}

					tokens[nTokenCount].pop_back();


					triangles.push_back({ verts[stoi(tokens[0]) - 1], verts[stoi(tokens[2]) - 1], verts[stoi(tokens[4]) - 1],
						texs[stoi(tokens[1]) - 1], texs[stoi(tokens[3]) - 1], texs[stoi(tokens[5]) - 1] });

				}

			}
		}
		return true;
	}
};

struct matrix4x4
{
	float matrix[4][4] = { 0 };
};

class Engine3D : public olcConsoleGameEngine
{
public:
	Engine3D()
	{
		m_sAppName = L"3D Demo";
	}

	bool OnUserCreate() override
	{
		pDepthBuffer = new float[ScreenWidth() * ScreenHeight()];

		sprTex1 = new olcSprite(L"test.spr"); //User change the argument to a custom spr file to render textures
		cube.loadObj("test.obj", true); //User change the argument to a custom obj file to render object

		projection_matrix = matrixProject(90.0f, (float)ScreenHeight() / (float)ScreenWidth(), 0.1f, 1000.0f);

		return true;
	}

	bool OnUserUpdate(float elapsedTime) override
	{
		if (GetKey(VK_UP).bHeld) camera.y += 8.0f * elapsedTime;
		if (GetKey(VK_DOWN).bHeld) camera.y -= 8.0f * elapsedTime;
		if (GetKey(VK_LEFT).bHeld) camera.x -= 8.0f * elapsedTime;
		if (GetKey(VK_RIGHT).bHeld) camera.x += 8.0f * elapsedTime;
		
		vec3d forward = vectorMul(look, 8.0f * elapsedTime);

		if (GetKey(L'W').bHeld) camera = vectorAdd(camera, forward);
		if (GetKey(L'S').bHeld) camera = vectorSub(camera, forward);
		if (GetKey(L'A').bHeld) yaw -= 2.0f * elapsedTime;
		if (GetKey(L'D').bHeld) yaw += 2.0f * elapsedTime;

		Fill(0, 0, ScreenWidth(), ScreenHeight(), PIXEL_SOLID, FG_BLACK);

		matrix4x4 rotateX, rotateZ;

		rotateX = matrixRotateX(theta);
		rotateZ = matrixRotateY(theta * 0.5);

		matrix4x4 translated;
		translated = matrixTranslate(0.0f, 0.0f, 16.0f);

		matrix4x4 world;
		world = matrixIdentity();
		world = matrixMultiply(rotateZ, rotateX);
		world = matrixMultiply(world, translated);

		vec3d up = { 0, 1, 0 };
		vec3d target = { 0, 0, 1 };
		matrix4x4 cameraRotationMatrix = matrixRotateY(yaw);
		look = matrixMultiplyVector(cameraRotationMatrix, target);
		target = vectorAdd(camera, look);

		matrix4x4 cameraMatrix = matrixPointAt(camera, target, up);
		matrix4x4 cameraView = matrixQuickInverse(cameraMatrix);

		vector<triangle> raster;

		for (auto tri : cube.triangles)
		{
			triangle projected, transformed, viewed;

			for (int i = 0; i < 3; i++)
			{
				transformed.p[i] = matrixMultiplyVector(world, tri.p[i]);
				transformed.t[i] = tri.t[i];
			}

			vec3d normal, line1, line2;
			line1 = vectorSub(transformed.p[1], transformed.p[0]);
			line2 = vectorSub(transformed.p[2], transformed.p[0]);
			normal = vectorCrossProduct(line1, line2);
			normal = vectorNormalise(normal);

			vec3d cameraRay = vectorSub(transformed.p[0], camera);

			if (calcDotProduct(normal, cameraRay) < 0.0f)
			{
				vec3d light = { 0.0f, 1.0f, -1.0f };
				light = vectorNormalise(light);
				float dot_product = max(0.1f, calcDotProduct(light, normal));

				CHAR_INFO c = GetColour(dot_product);
				transformed.col = c.Attributes;
				transformed.sym = c.Char.UnicodeChar;

				for (int i = 0; i < 3; i++)
				{
					viewed.p[i] = matrixMultiplyVector(cameraView, transformed.p[i]);
					viewed.t[i] = transformed.t[i];
				}
				viewed.sym = transformed.sym;
				viewed.col = transformed.col;

				int clippedTriangles = 0;
				triangle clipped[2];
				clippedTriangles = triangleClip({ 0.0f, 0.0f, 0.1f }, { 0.0f, 0.0f, 1.0f }, viewed, clipped[0], clipped[1]);

				for (int n = 0; n < clippedTriangles; n++)
				{
					for (int i = 0; i < 3; i++)
					{
						projected.p[i] = matrixMultiplyVector(projection_matrix, clipped[n].p[i]);
						projected.t[i] = clipped[n].t[i];
					}
					projected.col = clipped[n].col;
					projected.sym = clipped[n].sym;

					for (int i = 0; i < 3; i++)
					{
						projected.t[i].u = projected.t[i].u / projected.p[i].w;
						projected.t[i].v = projected.t[i].v / projected.p[i].w;
						projected.t[i].w = 1.0f / projected.p[i].w;
					}

					for (int i = 0; i < 3; i++)
					{
						projected.p[i] = vectorDiv(projected.p[i], projected.p[i].w);
						projected.p[i].x *= -1.0f;
						projected.p[i].y *= -1.0f;
					}

					vec3d offset = { 1, 1, 0 };
					for (int i = 0; i < 3; i++)
					{
						projected.p[i] = vectorAdd(projected.p[i], offset);
						projected.p[i].x *= 0.5f * (float)ScreenWidth();
						projected.p[i].y *= 0.5f * (float)ScreenHeight();
					}

					raster.push_back(projected);
				}
			}
		}

		sort(raster.begin(), raster.end(), [](triangle& t1, triangle& t2)
											{
												float z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z) / 3.0f;
												float z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z) / 3.0f;
												return z1 > z2;
											});

		Fill(0, 0, ScreenWidth(), ScreenHeight(), PIXEL_SOLID, FG_CYAN);

		for (int i = 0; i < ScreenWidth() * ScreenHeight(); i++)
		{
			pDepthBuffer[i] = 0.0f;
		}

		for (auto& toRaster : raster)
		{
			triangle clipped[2];
			list<triangle> listTriangles;
			listTriangles.push_back(toRaster);
			int newTriangles = 1;

			for (int p = 0; p < 4; p++)
			{
				int toAdd = 0;
				while (newTriangles > 0)
				{
					triangle test = listTriangles.front();
					listTriangles.pop_front();
					newTriangles--;

					switch (p)
					{
					case 0:	toAdd = triangleClip({ 0.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f }, test, clipped[0], clipped[1]); break;
					case 1:	toAdd = triangleClip({ 0.0f, (float)ScreenHeight() - 1, 0.0f }, { 0.0f, -1.0f, 0.0f }, test, clipped[0], clipped[1]); break;
					case 2:	toAdd = triangleClip({ 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;
					case 3:	toAdd = triangleClip({ (float)ScreenWidth() - 1, 0.0f, 0.0f }, { -1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;
					}

					for (int w = 0; w < toAdd; w++)
						listTriangles.push_back(clipped[w]);
				}
				newTriangles = listTriangles.size();
			}

			for (auto& t : listTriangles)
			{
				drawTexturedTriangle(t.p[0].x, t.p[0].y, t.t[0].u, t.t[0].v, t.t[0].w,
					             t.p[1].x, t.p[1].y, t.t[1].u, t.t[1].v, t.t[1].w,
					             t.p[2].x, t.p[2].y, t.t[2].u, t.t[2].v, t.t[2].w, sprTex1);
			}
		}

		return true;
	}

	void drawTexturedTriangle(int x1, int y1, float u1, float v1, float w1, int x2, int y2, float u2, float v2, float w2, int x3, int y3, float u3, float v3, float w3, olcSprite *tex)
	{
		if (y2 < y1)
		{
			swap(y1, y2);
			swap(x1, x2);
			swap(u1, u2);
			swap(v1, v2);
			swap(w1, w2);
		}

		if (y3 < y1)
		{
			swap(y1, y3);
			swap(x1, x3);
			swap(u1, u3);
			swap(v1, v3);
			swap(w1, w3);
		}

		if (y3 < y2)
		{
			swap(y2, y3);
			swap(x2, x3);
			swap(u2, u3);
			swap(v2, v3);
			swap(w2, w3);
		}

		int dy1 = y2 - y1;
		int dx1 = x2 - x1;
		float dv1 = v2 - v1;
		float du1 = u2 - u1;
		float dw1 = w2 - w1;

		int dy2 = y3 - y1;
		int dx2 = x3 - x1;
		float dv2 = v3 - v1;
		float du2 = u3 - u1;
		float dw2 = w3 - w1;

		float tex_u, tex_v, tex_w;

		float dax_step = 0, dbx_step = 0,
			du1_step = 0, dv1_step = 0,
			du2_step = 0, dv2_step = 0,
			dw1_step = 0, dw2_step = 0;

		if (dy1) dax_step = dx1 / (float)abs(dy1);
		if (dy2) dbx_step = dx2 / (float)abs(dy2);

		if (dy1) du1_step = du1 / (float)abs(dy1);
		if (dy1) dv1_step = dv1 / (float)abs(dy1);
		if (dy1) dw1_step = dw1 / (float)abs(dy1);

		if (dy2) du2_step = du2 / (float)abs(dy2);
		if (dy2) dv2_step = dv2 / (float)abs(dy2);
		if (dy2) dw2_step = dw2 / (float)abs(dy2);

		if (dy1)
		{
			for (int i = y1; i <= y2; i++)
			{
				int ax = x1 + (float)(i - y1) * dax_step;
				int bx = x1 + (float)(i - y1) * dbx_step;

				float tex_su = u1 + (float)(i - y1) * du1_step;
				float tex_sv = v1 + (float)(i - y1) * dv1_step;
				float tex_sw = w1 + (float)(i - y1) * dw1_step;

				float tex_eu = u1 + (float)(i - y1) * du2_step;
				float tex_ev = v1 + (float)(i - y1) * dv2_step;
				float tex_ew = w1 + (float)(i - y1) * dw2_step;

				if (ax > bx)
				{
					swap(ax, bx);
					swap(tex_su, tex_eu);
					swap(tex_sv, tex_ev);
					swap(tex_sw, tex_ew);
				}

				tex_u = tex_su;
				tex_v = tex_sv;
				tex_w = tex_sw;

				float tstep = 1.0f / ((float)(bx - ax));
				float t = 0.0f;

				for (int j = ax; j < bx; j++)
				{
					tex_u = (1.0f - t) * tex_su + t * tex_eu;
					tex_v = (1.0f - t) * tex_sv + t * tex_ev;
					tex_w = (1.0f - t) * tex_sw + t * tex_ew;
					if (tex_w > pDepthBuffer[i * ScreenWidth() + j])
					{
						Draw(j, i, tex->SampleGlyph(tex_u / tex_w, tex_v / tex_w), tex->SampleColour(tex_u / tex_w, tex_v / tex_w));
						pDepthBuffer[i * ScreenWidth() + j] = tex_w;
					}
					t += tstep;
				}
			}
		}

		dy1 = y3 - y2;
		dx1 = x3 - x2;
		dv1 = v3 - v2;
		du1 = u3 - u2;
		dw1 = w3 - w2;

		if (dy1) dax_step = dx1 / (float)abs(dy1);
		if (dy2) dbx_step = dx2 / (float)abs(dy2);

		du1_step = 0, dv1_step = 0;
		if (dy1) du1_step = du1 / (float)abs(dy1);
		if (dy1) dv1_step = dv1 / (float)abs(dy1);
		if (dy1) dw1_step = dw1 / (float)abs(dy1);

		if (dy1)
		{
			for (int i = y2; i <= y3; i++)
			{
				int ax = x2 + (float)(i - y2) * dax_step;
				int bx = x1 + (float)(i - y1) * dbx_step;

				float tex_su = u2 + (float)(i - y2) * du1_step;
				float tex_sv = v2 + (float)(i - y2) * dv1_step;
				float tex_sw = w2 + (float)(i - y2) * dw1_step;

				float tex_eu = u1 + (float)(i - y1) * du2_step;
				float tex_ev = v1 + (float)(i - y1) * dv2_step;
				float tex_ew = w1 + (float)(i - y1) * dw2_step;

				if (ax > bx)
				{
					swap(ax, bx);
					swap(tex_su, tex_eu);
					swap(tex_sv, tex_ev);
					swap(tex_sw, tex_ew);
				}

				tex_u = tex_su;
				tex_v = tex_sv;
				tex_w = tex_sw;

				float tstep = 1.0f / ((float)(bx - ax));
				float t = 0.0f;

				for (int j = ax; j < bx; j++)
				{
					tex_u = (1.0f - t) * tex_su + t * tex_eu;
					tex_v = (1.0f - t) * tex_sv + t * tex_ev;
					tex_w = (1.0f - t) * tex_sw + t * tex_ew;

					if (tex_w > pDepthBuffer[i * ScreenWidth() + j])
					{
						Draw(j, i, tex->SampleGlyph(tex_u / tex_w, tex_v / tex_w), tex->SampleColour(tex_u / tex_w, tex_v / tex_w));
						pDepthBuffer[i * ScreenWidth() + j] = tex_w;
					}
					t += tstep;
				}
			}
		}
	}

private:
	mesh cube;
	matrix4x4 projection_matrix;
	float theta;
	
	vec3d camera;
	vec3d look;

	float yaw;

	olcSprite* sprTex1;
	float* pDepthBuffer = nullptr;

	vec3d matrixMultiplyVector(matrix4x4& m, vec3d& i)
	{
		vec3d v;
		v.x = i.x * m.matrix[0][0] + i.y * m.matrix[1][0] + i.z * m.matrix[2][0] + i.w * m.matrix[3][0];
		v.y = i.x * m.matrix[0][1] + i.y * m.matrix[1][1] + i.z * m.matrix[2][1] + i.w * m.matrix[3][1];
		v.z = i.x * m.matrix[0][2] + i.y * m.matrix[1][2] + i.z * m.matrix[2][2] + i.w * m.matrix[3][2];
		v.w = i.x * m.matrix[0][3] + i.y * m.matrix[1][3] + i.z * m.matrix[2][3] + i.w * m.matrix[3][3];
		return v;
	}

	matrix4x4 matrixIdentity()
	{
		matrix4x4 matrix;
		matrix.matrix[0][0] = 1.0f;
		matrix.matrix[1][1] = 1.0f;
		matrix.matrix[2][2] = 1.0f;
		matrix.matrix[3][3] = 1.0f;
		return matrix;
	}

	matrix4x4 matrixRotateX(float fAngleRad)
	{
		matrix4x4 matrix;
		matrix.matrix[0][0] = 1.0f;
		matrix.matrix[1][1] = cosf(fAngleRad);
		matrix.matrix[1][2] = sinf(fAngleRad);
		matrix.matrix[2][1] = -sinf(fAngleRad);
		matrix.matrix[2][2] = cosf(fAngleRad);
		matrix.matrix[3][3] = 1.0f;
		return matrix;
	}

	matrix4x4 matrixRotateY(float fAngleRad)
	{
		matrix4x4 matrix;
		matrix.matrix[0][0] = cosf(fAngleRad);
		matrix.matrix[0][2] = sinf(fAngleRad);
		matrix.matrix[2][0] = -sinf(fAngleRad);
		matrix.matrix[1][1] = 1.0f;
		matrix.matrix[2][2] = cosf(fAngleRad);
		matrix.matrix[3][3] = 1.0f;
		return matrix;
	}

	matrix4x4 matrixRotateZ(float fAngleRad)
	{
		matrix4x4 matrix;
		matrix.matrix[0][0] = cosf(fAngleRad);
		matrix.matrix[0][1] = sinf(fAngleRad);
		matrix.matrix[1][0] = -sinf(fAngleRad);
		matrix.matrix[1][1] = cosf(fAngleRad);
		matrix.matrix[2][2] = 1.0f;
		matrix.matrix[3][3] = 1.0f;
		return matrix;
	}

	matrix4x4 matrixTranslate(float x, float y, float z)
	{
		matrix4x4 matrix;
		matrix.matrix[0][0] = 1.0f;
		matrix.matrix[1][1] = 1.0f;
		matrix.matrix[2][2] = 1.0f;
		matrix.matrix[3][3] = 1.0f;
		matrix.matrix[3][0] = x;
		matrix.matrix[3][1] = y;
		matrix.matrix[3][2] = z;
		return matrix;
	}

	matrix4x4 matrixProject(float fFovDegrees, float fAspectRatio, float fNear, float fFar)
	{
		float fFovRad = 1.0f / tanf(fFovDegrees * 0.5f / 180.0f * 3.14159f);
		matrix4x4 matrix;
		matrix.matrix[0][0] = fAspectRatio * fFovRad;
		matrix.matrix[1][1] = fFovRad;
		matrix.matrix[2][2] = fFar / (fFar - fNear);
		matrix.matrix[3][2] = (-fFar * fNear) / (fFar - fNear);
		matrix.matrix[2][3] = 1.0f;
		matrix.matrix[3][3] = 0.0f;
		return matrix;
	}

	matrix4x4 matrixMultiply(matrix4x4& m1, matrix4x4& m2)
	{
		matrix4x4 matrix;
		for (int c = 0; c < 4; c++)
			for (int r = 0; r < 4; r++)
				matrix.matrix[r][c] = m1.matrix[r][0] * m2.matrix[0][c] + m1.matrix[r][1] * m2.matrix[1][c] + m1.matrix[r][2] * m2.matrix[2][c] + m1.matrix[r][3] * m2.matrix[3][c];
		return matrix;
	}

	matrix4x4 matrixPointAt(vec3d& pos, vec3d& target, vec3d& up)
	{
		vec3d newForward = vectorSub(target, pos);
		newForward = vectorNormalise(newForward);

		vec3d a = vectorMul(newForward, calcDotProduct(up, newForward));
		vec3d newUp = vectorSub(up, a);
		newUp = vectorNormalise(newUp);

		vec3d newRight = vectorCrossProduct(newUp, newForward);
	
		matrix4x4 matrix;
		matrix.matrix[0][0] = newRight.x;	matrix.matrix[0][1] = newRight.y;	matrix.matrix[0][2] = newRight.z;	matrix.matrix[0][3] = 0.0f;
		matrix.matrix[1][0] = newUp.x;		matrix.matrix[1][1] = newUp.y;		matrix.matrix[1][2] = newUp.z;		matrix.matrix[1][3] = 0.0f;
		matrix.matrix[2][0] = newForward.x;	matrix.matrix[2][1] = newForward.y;	matrix.matrix[2][2] = newForward.z;	matrix.matrix[2][3] = 0.0f;
		matrix.matrix[3][0] = pos.x;		matrix.matrix[3][1] = pos.y;		matrix.matrix[3][2] = pos.z;		matrix.matrix[3][3] = 1.0f;
		return matrix;

	}

	matrix4x4 matrixQuickInverse(matrix4x4& m)
	{
		matrix4x4 matrix;
		matrix.matrix[0][0] = m.matrix[0][0]; matrix.matrix[0][1] = m.matrix[1][0]; matrix.matrix[0][2] = m.matrix[2][0]; matrix.matrix[0][3] = 0.0f;
		matrix.matrix[1][0] = m.matrix[0][1]; matrix.matrix[1][1] = m.matrix[1][1]; matrix.matrix[1][2] = m.matrix[2][1]; matrix.matrix[1][3] = 0.0f;
		matrix.matrix[2][0] = m.matrix[0][2]; matrix.matrix[2][1] = m.matrix[1][2]; matrix.matrix[2][2] = m.matrix[2][2]; matrix.matrix[2][3] = 0.0f;
		matrix.matrix[3][0] = -(m.matrix[3][0] * matrix.matrix[0][0] + m.matrix[3][1] * matrix.matrix[1][0] + m.matrix[3][2] * matrix.matrix[2][0]);
		matrix.matrix[3][1] = -(m.matrix[3][0] * matrix.matrix[0][1] + m.matrix[3][1] * matrix.matrix[1][1] + m.matrix[3][2] * matrix.matrix[2][1]);
		matrix.matrix[3][2] = -(m.matrix[3][0] * matrix.matrix[0][2] + m.matrix[3][1] * matrix.matrix[1][2] + m.matrix[3][2] * matrix.matrix[2][2]);
		matrix.matrix[3][3] = 1.0f;
		return matrix;
	}


	vec3d vectorAdd(vec3d& v1, vec3d& v2)
	{
		return { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
	}

	vec3d vectorSub(vec3d& v1, vec3d& v2)
	{
		return { v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };
	}

	vec3d vectorMul(vec3d& v1, float k)
	{
		return { v1.x * k, v1.y * k, v1.z * k };
	}

	vec3d vectorDiv(vec3d& v1, float k)
	{
		return { v1.x / k, v1.y / k, v1.z / k };
	}

	float calcDotProduct(vec3d& v1, vec3d& v2)
	{
		return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
	}

	float vectorLength(vec3d& v)
	{
		return sqrtf(calcDotProduct(v, v));
	}

	vec3d vectorNormalise(vec3d& v)
	{
		float l = vectorLength(v);
		return { v.x / l, v.y / l, v.z / l };
	}

	vec3d vectorCrossProduct(vec3d& v1, vec3d& v2)
	{
		vec3d v;
		v.x = v1.y * v2.z - v1.z * v2.y;
		v.y = v1.z * v2.x - v1.x * v2.z;
		v.z = v1.x * v2.y - v1.y * v2.x;
		return v;
	}

	vec3d vectorIntersectPlane(vec3d& plane_p, vec3d& plane_n, vec3d& lineStart, vec3d& lineEnd, float &t)
	{
		plane_n = vectorNormalise(plane_n);
		float plane_d = -calcDotProduct(plane_n, plane_p);
		float ad = calcDotProduct(lineStart, plane_n);
		float bd = calcDotProduct(lineEnd, plane_n);
		t = (-plane_d - ad) / (bd - ad);
		vec3d lineStartToEnd = vectorSub(lineEnd, lineStart);
		vec3d lineToIntersect = vectorMul(lineStartToEnd, t);
		return vectorAdd(lineStart, lineToIntersect);
	}

	int triangleClip(vec3d plane_p, vec3d plane_n, triangle& in_tri, triangle& out_tri1, triangle& out_tri2)
	{
		plane_n = vectorNormalise(plane_n);

		auto dist = [&](vec3d& p)
		{
			vec3d n = vectorNormalise(p);
			return (plane_n.x * p.x + plane_n.y * p.y + plane_n.z * p.z - calcDotProduct(plane_n, plane_p));
		};

		vec3d* inside_points[3];  int nInsidePointCount = 0;
		vec3d* outside_points[3]; int nOutsidePointCount = 0;

		vec2d* inside_tex[3]; int nInsideTexCount = 0;
		vec2d* outside_tex[3]; int nOutsideTexCount = 0;

		float d0 = dist(in_tri.p[0]);
		float d1 = dist(in_tri.p[1]);
		float d2 = dist(in_tri.p[2]);

		if (d0 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[0]; inside_tex[nInsideTexCount++] = &in_tri.t[0]; }
		else { outside_points[nOutsidePointCount++] = &in_tri.p[0]; outside_tex[nOutsideTexCount++] = &in_tri.t[0]; }

		if (d1 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[1]; inside_tex[nInsideTexCount++] = &in_tri.t[1]; }
		else { outside_points[nOutsidePointCount++] = &in_tri.p[1];  outside_tex[nOutsideTexCount++] = &in_tri.t[1]; }

		if (d2 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[2]; inside_tex[nInsideTexCount++] = &in_tri.t[2]; }
		else { outside_points[nOutsidePointCount++] = &in_tri.p[2];  outside_tex[nOutsideTexCount++] = &in_tri.t[2]; }

		if (nInsidePointCount == 0) { return 0; }

		if (nInsidePointCount == 3) { out_tri1 = in_tri; return 1; }

		if (nInsidePointCount == 1 && nOutsidePointCount == 2)
		{
			out_tri1.col = in_tri.col;
			out_tri1.sym = in_tri.sym;

			out_tri1.p[0] = *inside_points[0];
			out_tri1.t[0] = *inside_tex[0];

			float t;
			out_tri1.p[1] = vectorIntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0], t);
			out_tri1.t[1].u = t * (outside_tex[0]->u - inside_tex[0]->u) + inside_tex[0]->u;
			out_tri1.t[1].v = t * (outside_tex[0]->v - inside_tex[0]->v) + inside_tex[0]->v;
			out_tri1.t[1].w = t * (outside_tex[0]->w - inside_tex[0]->w) + inside_tex[0]->w;

			out_tri1.p[2] = vectorIntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[1], t);
			out_tri1.t[2].u = t * (outside_tex[1]->u - inside_tex[0]->u) + inside_tex[0]->u;
			out_tri1.t[2].v = t * (outside_tex[1]->v - inside_tex[0]->v) + inside_tex[0]->v;
			out_tri1.t[2].w = t * (outside_tex[1]->w - inside_tex[0]->w) + inside_tex[0]->w;

			return 1;
		}

		if (nInsidePointCount == 2 && nOutsidePointCount == 1)
		{
			out_tri1.col = in_tri.col;
			out_tri1.sym = in_tri.sym;

			out_tri2.col = in_tri.col;
			out_tri2.sym = in_tri.sym;

			out_tri1.p[0] = *inside_points[0];
			out_tri1.p[1] = *inside_points[1];
			out_tri1.t[0] = *inside_tex[0];
			out_tri1.t[1] = *inside_tex[1];

			float t;
			out_tri1.p[2] = vectorIntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0], t);
			out_tri1.t[2].u = t * (outside_tex[0]->u - inside_tex[0]->u) + inside_tex[0]->u;
			out_tri1.t[2].v = t * (outside_tex[0]->v - inside_tex[0]->v) + inside_tex[0]->v;
			out_tri1.t[2].w = t * (outside_tex[0]->w - inside_tex[0]->w) + inside_tex[0]->w;

			out_tri2.p[0] = *inside_points[1];
			out_tri2.t[0] = *inside_tex[1];
			out_tri2.p[1] = out_tri1.p[2];
			out_tri2.t[1] = out_tri1.t[2];
			out_tri2.p[2] = vectorIntersectPlane(plane_p, plane_n, *inside_points[1], *outside_points[0], t);
			out_tri2.t[2].u = t * (outside_tex[0]->u - inside_tex[1]->u) + inside_tex[1]->u;
			out_tri2.t[2].v = t * (outside_tex[0]->v - inside_tex[1]->v) + inside_tex[1]->v;
			out_tri2.t[2].w = t * (outside_tex[0]->w - inside_tex[1]->w) + inside_tex[1]->w;

			return 2;
		}
	}

	CHAR_INFO GetColour(float lum)
	{
		short bg_col, fg_col;
		wchar_t sym;
		int pixel_bw = (int)(13.0f * lum);

		switch (pixel_bw)
		{
		case 0: bg_col = BG_BLACK; fg_col = FG_BLACK; sym = PIXEL_SOLID; break;

		case 1: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_QUARTER; break;
		case 2: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_HALF; break;
		case 3: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_THREEQUARTERS; break;
		case 4: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_SOLID; break;

		case 5: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_QUARTER; break;
		case 6: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_HALF; break;
		case 7: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_THREEQUARTERS; break;
		case 8: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_SOLID; break;

		case 9:  bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_QUARTER; break;
		case 10: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_HALF; break;
		case 11: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_THREEQUARTERS; break;
		case 12: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_SOLID; break;
		default:
			bg_col = BG_BLACK; fg_col = FG_BLACK; sym = PIXEL_SOLID;
		}

		CHAR_INFO c;
		c.Attributes = bg_col | fg_col;
		c.Char.UnicodeChar = sym;
		return c;
	}
};


int main(int argc, char* argv[])
{
	Engine3D demo;

	if (demo.ConstructConsole(256, 240, 4, 4))
	{
		demo.Start();
	}

	return 0;
}