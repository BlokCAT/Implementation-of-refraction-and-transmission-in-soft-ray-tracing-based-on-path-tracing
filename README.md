# 一种在基于路径追踪（Path_Tracing）的软光追渲染器中实现透射和折射的有效实现方案



软光追渲染器Demo的链接：

[BlokCAT/Ray-Tracing-software-renderer: 使用纯C++和文件输出流实现的一个基于蒙德卡罗积分路径追踪的软光线追踪渲染器，结构清晰，体量中等，比较适合学习参考 (github.com)](https://github.com/BlokCAT/Ray-Tracing-software-renderer)

**建议先了解实现软光追原理过程再看，这里只是介绍实现折射的原理，完整代码在上面demo！！**

## 展示

上图：

![image](https://github.com/BlokCAT/Implementation-of-refraction-and-transmission-in-soft-ray-tracing-based-on-path-tracing/blob/main/1.png?raw=true)

其中有聚焦效果，说明实现没有问题：

![image-20240921120641537](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20240921120641537.png)

## 原理详细介绍

主要需要看下面的图

<img src="README.assets/image-20240921132650615.png" alt="image-20240921132650615" style="zoom:50%;" />

最终的目的是计算出`Color(light_in)`的颜色,而这个颜色来源是反射颜色`Color(reflect)`和折射颜色`Color(refract_light)`的加权和，只需要用菲尼尔效应计算出各自的占比就行了，

​	1.反射光使用镜面反射的方式直接计算得出（项目里有详细代码）

​	2.折射光refract进去后需要重新找交点， 在出去的点再次计算出折射方向，形成一个新光线继续递归计算颜色，就得到了`Color(refract_light)`的值（注意向量都是由着色点为起点的）

**所以材质类里面的一些函数如下：**

菲尼尔:

```c
	void fresnel(const Vector3f &in, const Vector3f &N, float &ior, float &kr) //kr是反射的比例
	{
		in.normalized();
		N.normalized();
		float Ro = ((1.0 - ior) / (1.0 + ior)) * ((1.0 - ior) / (1.0 + ior));
		float costheta = clamp( 0.0f , 1.0 ,  dotProduct(N, in));
		float t = (1.0 - costheta) * (1.0 - costheta);
		kr = Ro + (1 - Ro) * t;
		return;
	}
```

计算折射方向（考虑了全反射和射出和射入的方向）

```c
	Vector3f refract(const Vector3f &II, const Vector3f &N, const float &ior) 
	{
		Vector3f I = II * -1;
		// 计算入射角的余弦值
		float cosi = clamp(-1, 1, dotProduct(I, N));
		// 确定折射率比和法线方向
		float etai = 1, etat = ior;
		Vector3f n = N;
		if (cosi < 0)cosi = -cosi; // 入射光线从外部进入物体内部
		else{std::swap(etai, etat); n = N * -1;}
		// 计算折射率比
		float eta = etai / etat;
		// 计算折射后的光线方向
		float k = 1 - eta * eta * (1 - cosi * cosi);
		if (k < 0) return Vector3f(0.0f);
		else{
			// 计算折射光线的方向
			return (I* eta + ( n * (eta * cosi - sqrtf(k)))).normalized();
		}
	}
```

`getBRDF`需要加一个新的，因为射入时需要表面的材质颜色和射出不需要，所以要分开

```c
Vector3f GetRefracBRDF(const Vector3f& wi, const Vector3f& wo, const Vector3f& N , int flag)
{
	switch (flag)
	{
	case 1:
	{
		float ans = dotProduct(wo, N);
		if (ans > 0.0f)
		{
			Vector3f res = 1.0 / ans;
			return res;
		}
		return Vector3f(0);
	}
	case 2:
	{
		float ans = dotProduct(wo, N);
		if (ans > 0.0f)
		{
			Vector3f res = 1.0 / ans;
			res = (res * Kd);
			return res;
		}
		return Vector3f(0);
		break;
	}
	default:
		break;
	}

}
```
PDF的话， 如果角度正确就返回1

```c
float pdf(const Vector3f &wi, const Vector3f &wo, const Vector3f &N)
{
	switch (this->mtype)
	{
	case MIRCO:
	case DIFFUSE:
	{
		float ans = dotProduct(wo, N);
		if (ans > 0.0f)
		{
			return 0.5f / M_PI;
		}
		return 0.0f;
		break;
	}
	case REFRACT:  // 这个是折射的pdf
	{
		float ans = dotProduct(wo, N);
		if (ans > 0.0)
		{
			return 1.0;
		}
		return 0.0f;

		break;
	}
	case REFLC:
	{
		float ans = dotProduct(wo, N);
		if (ans > 0.0)
		{
			return 1.0;
		}
		return 0.0f;
		break;
	}
	default:
		break;
	}
}

```

其中操作这些函数的`Scene.cpp`类

```c

//Path Tracing
Vector3f Scene::PathTracing( Ray &ray, int depth)
{

	Vector3f L_dir(0.0), L_indir(0.0);
	//if (depth >= 1)return L_dir; 

	HitPoint hit_to_scene;
	hit_to_scene.happened = false;
	Scene::FindHit(ray , hit_to_scene);

	
	if (!hit_to_scene.happened) return Vector3f(0);
	if (hit_to_scene.m->islight) return hit_to_scene.m->lightIntensity;

	//取出所有光线在场景的第一个交点的信息
	Vector3f hit_pos = hit_to_scene.hitcoord;
	Vector3f N = hit_to_scene.hitN.normalized();
	Material * mat = hit_to_scene.m;
	Vector3f wi = (ray.dir * -1).normalized();

	switch (mat->mtype)
	{
	case REFLC:
	{
		// 详细看项目代码
	}
	case MIRCO:
	case DIFFUSE:
	{
		// 详细看项目代码			
	}
	case REFRACT:  //折射代码
	{
		float gs = RandomFloat();
		if( gs < RussianRoulette)
		{
			Vector3f Color_rflcet(0), Color_rfract(0);
			float K = 0.0;
			mat->fresnel(wi, N, mat->ior, K);

			//计算折射颜色
			Vector3f futureDir1 = mat->refract(wi, N , mat->ior);
			Ray newRay1(hit_pos - (N * 0.001), futureDir1);
			HitPoint test1;
			test1.happened = false;
			FindHit(newRay1, test1);
	
			if (test1.happened)
			{
				//取出所有光线在出去的点的信息
				Vector3f hitpos = test1.hitcoord;
				Vector3f nn = test1.hitN.normalized();
				Material* mt = test1.m;
				Vector3f new_wi = (newRay1.dir * -1).normalized();
				Vector3f futureOutDir = mt->refract(new_wi, nn, mt->ior);

				if (!(!futureOutDir.x && !futureOutDir.y && !futureOutDir.z)) {
					//折射后射出去的光线
	
					Ray outRay(hitpos + (nn * 0.001), futureOutDir);
					HitPoint Obj_rfract_hit;
					FindHit(outRay, Obj_rfract_hit);

					if (Obj_rfract_hit.happened)
					{
						Vector3f Brdf_ = mt->GetRefracBRDF(new_wi, futureOutDir, nn , 1);
						float costheta4 = dotProduct(futureOutDir, nn);
						float Pdf_ = mt->pdf(new_wi, futureOutDir, nn);
						if (Pdf_ > 0.0001)
							Color_rfract = PathTracing(outRay, depth + 1) * Brdf_ * costheta4 / Pdf_ / RussianRoulette;
					}
				}
			
			}
			
			//计算折射体表面反射的颜色
			Vector3f futureDir2 = mat->GetFutureDir(wi, N);
			Ray newRay2(hit_pos, futureDir2);
			HitPoint test2;
			FindHit(newRay2, test2);
			if (test2.happened)
			{
				Vector3f brdf_ = mat->GetRefracBRDF(wi, futureDir2, N , 2);
				float costheta5 = dotProduct(futureDir2, N);
				float pdf_ = mat->pdf(wi, futureDir2, N);
				if (pdf_ > 0.0001)
					Color_rflcet = PathTracing(newRay2, depth + 1) * brdf_ * costheta5 / pdf_ / RussianRoulette;
			}
			L_indir = (Color_rfract) * (1 - K) + (Color_rflcet)*K;
		}

		break;

	}
	default:
		break;
	}
	
	Vector3f res = L_dir + L_indir;
	res.x = clamp(0.0, 14, res.x);
	res.y = clamp(0.0, 14, res.y);
	res.z = clamp(0.0, 14, res.z);
	return res;
}
	
```



