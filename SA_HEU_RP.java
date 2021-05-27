package TCybernetics;

import java.io.File;
import java.io.FileInputStream;
import java.util.Date;
import jxl.Workbook;
import jxl.write.WritableSheet;
import jxl.write.WritableWorkbook;

/*
 * 3个benchmark函数：1)transmit ;2)send ;3)processTupleArrival ;4)executeTuple ;5)checkCloudletCompletion ;6)getRsultantTuples
 * (函数序号，输入变量维度，存在路径数目，可行路径数目)
 * transmit(1,3,2,1) 100%; send(2,2,9,5)66%; processEvent(3,7,9,6)100%; executeTuple(4,7,5,3)100%; 
 * checkCloudletCompletion(5,5,6,3)100%; getResultantTuple(6,8,7,4)100%
 * getPhysicalTopology(7,110,7,6);  gimp_rgb_to_hsv(8,3,6,5);  gimp_hsv_to_rgb(9,3,15,3)
 * check_ISBN(10,12,7,5);  check_ISSN(11,10,7,5);  gimp_hwb_to_rgb(12,3,15,3)
 * TIFF_GetSourceSamples(13,2,8,7)
 */

public class SA_HEU_RP {
	private static final int RUN = 30;                       //运行次数
	private static final int pop_num = 50;                   //种群大小
	
	/*适应值计算参数*/
	private static final int K = 10;               
	private static final double alpha = 0.001;    
	
	/*DE算法参数*/
	private static final double Pc = 0.2;
	private static final double F = 0.5;
	
	/*测试函数信息*/
	private static int fun_num;
	private static int R;
	private static int PATHNUM;
	private static int NODENUM;
	private static final int BRANCH = 15;
	private static boolean[][] visit;
	
	/* infection是一个n*m的矩阵C
	 * 矩阵C负责记录测试用例的第i个维度，对节点j的影响次数*/
	private static int[][]infection;
	/* record负责记录种群每个个体（测试用例）对应的节点覆盖情况*/
	private static int[] record = new int[pop_num];
	/*PATH记录不同路径的分支节点走向*/
	private static String[] PATH;
	/*搜索步长*/
	private static final int step_length1 = 2;
	private static final int step_length2 = 4;
	
	/*开始结束时间、运行时间、覆盖率、循环次数、测试用例数目*/
    static double start;                                             
    static double finish;
    static double[] runtime;
    static double[] coverage;
    static int[] case_num;
    static int[] obj = new int[1];
    
    /*输入变量参数上下界*/
	private static int[] lb;
	private static int[] ub;
	
	private static final int MCN = 300000;                  //最大迭代次数
	private static int col;
	
	public static void main(String[] args){
		
		for(fun_num = 13;fun_num<= 13;fun_num++)
		{		
			System.out.println("FUNCTION = " + fun_num);

			setFunctionParameters();// 函数基本参数设置

			lb = new int[R];
			ub = new int[R];
			visit = new boolean[NODENUM][BRANCH];
			PATH = new String[PATHNUM];
			infection = new int[R][NODENUM];
			runtime = new double[RUN];
			coverage = new double[RUN];
			case_num = new int[RUN];
			setDomainAndEncoding(); // 初始化函数定义域与路径编码

			for (int run = 0; run < RUN; run++) 
			{
				/* 初始化visit */
				for (int i = 0; i < NODENUM; i++)
					for (int j = 0; j < BRANCH; j++)
						visit[i][j] = false;
				/* 初始化infection */
				for (int i = 0; i < R; i++)
					for (int j = 0; j < NODENUM; j++)
						infection[i][j] = 0;

				/* 路径覆盖状态、状态初始化 */
				int[][] x = new int[pop_num][R];
				int[][] v = new int[pop_num][R];
				double[] fitness_x = new double[pop_num];
				double[] fitness_v = new double[pop_num];
				boolean[] status = new boolean[PATHNUM]; // to mark whether the path
															// has been covered.
				int[][] solution = new int[PATHNUM][R];
				int path;
				obj[0] = 0;

				Date mydate = new Date();
				start = mydate.getTime();

				/* 初始化种群 */
				for (int i = 0; i < pop_num; i++) {
					for (int j = 0; j < R; j++) {
						x[i][j] = (int) (lb[j] + Math.random() * (ub[j] - lb[j]));
					}

					if (obj[0] == PATHNUM)
						break;

					path = pathnum(x[i], fun_num); // 获取测试用例覆盖的路径
					record[i] = path; // 记录测试用例覆盖的路径编号

					if (!status[path]) {
						for (int j = 0; j < R; j++)
							solution[path][j] = x[i][j]; // 记录路径第一个
						status[path] = true; // 标记路径Path是否已找到覆盖它的用例
						obj[0]++; // 已覆盖的路径数
						nodeiscoverage(x[i], fun_num); // 标记已被覆盖的分支
					}
					case_num[run] = case_num[run] + 1; // 评估次数自增
				}
				
				/* Scatter Search搜索个体，但不更新关系矩阵 */

				ScatterSearch(x, fitness_x, solution, status, run, obj);
//				Homo_Euclidean_Update(x, fitness_x, solution, status, run, obj);
				
				System.out.println(case_num[run]);
				
				while (case_num[run] <= MCN && obj[0] < PATHNUM) 
				{
//					DE_search(x, v, fitness_x, fitness_v, solution, status, run, obj);
					
					if(obj[0] == PATHNUM)
						break;
					
					RelationshipProgramming(x, fitness_x, solution, status, run, obj);
				}
				
				Date mydate2 = new Date();
				finish = mydate2.getTime();
				runtime[run] = finish - start;
				System.out.println();
				System.out.println("运行时间=" + runtime[run] + "ms"); // 输出运行时间

				coverage[run] = obj[0] * 100.0 / PATHNUM;
				System.out.println("路径覆盖率=" + coverage[run] + "%");
				System.out.println("最优解为：");

				for (int k = 0; k < PATHNUM; k++) // 输出路径覆盖情况：覆盖路径的测试用例以及未覆盖的路径
				{
					if (status[k]) {
						System.out.print("path" + k + ":");
						for (int j = 0; j < R; j++)
							System.out.print((int) solution[k][j] + " ");
						System.out.println();
					} else
						System.out.println("path" + k + "没被覆盖.");
				}
				System.out.println("case_num[" + run + "] = " + case_num[run]);

				/* 输出关系矩阵 */
				System.out.println("\n" + "infection矩阵:");
				for (int r = 0; r < R; r++) {
					for (int n = 0; n < NODENUM; n++)
						System.out.print(infection[r][n] + "  ");
					System.out.println();
				}
			}

			double time_sum = 0, time_average, coverage_sum = 0, coverage_average, case_average;
			int case_sum = 0;
			for (int run = 0; run < RUN; run++) {
				time_sum = time_sum + runtime[run];
				coverage_sum = coverage_sum + coverage[run];
				case_sum = case_sum + case_num[run];
			}
			time_average = time_sum / RUN;
			coverage_average = coverage_sum / RUN;
			case_average = case_sum / RUN;

			System.out.println("time_sum = " + time_sum + "ms");
			System.out.println("time_average = " + time_average + "ms");
			System.out.println("case_sum = " + case_sum);
			System.out.println("case_average = " + case_average);
			System.out.println("coverage_sum = " + coverage_sum + "%");
			System.out.println("coverage_average = " + coverage_average + "%");

			try // 将数据导出到Excel文档中
			{
				File file = new java.io.File("G:/TCyber/temp.xls");
//				File file = new java.io.File("G:/TCyber/MIA.xls");
//				File file = new java.io.File("G:/TCyber/SA-HEU-RP.xls");
//				File file = new java.io.File("G:/TCyber/SS.xls");
//				FileInputStream in = new FileInputStream(file);
				Workbook book = Workbook.getWorkbook(file);
				WritableWorkbook wbook = Workbook.createWorkbook(file, book);
				WritableSheet sheet = wbook.getSheet(0); // 写入数据sheet

				jxl.write.Number ID = new jxl.write.Number(col, 0, fun_num);
				sheet.addCell(ID);
				
				for (int run = 0; run < RUN; run++) {
					int q = run;
					jxl.write.Number number = new jxl.write.Number(col, q+1,case_num[run]); 
					jxl.write.Number number2 = new jxl.write.Number(col, q+RUN+10,coverage[run]); 
					jxl.write.Number number3 = new jxl.write.Number(col, q+15+RUN*2, runtime[run]);
					sheet.addCell(number); 
					sheet.addCell(number2);
					sheet.addCell(number3);
				}

				double case_ave = getAverage(case_num, RUN);
				jxl.write.Number number1 = new jxl.write.Number(col, RUN + 4, case_ave);
				sheet.addCell(number1);
				double case_std = getStandardDevition(case_num, RUN);
				jxl.write.Number number2 = new jxl.write.Number(col, RUN + 5, case_std);
				sheet.addCell(number2);
				jxl.write.Number number3 = new jxl.write.Number(col, RUN + 6, coverage_average);
				sheet.addCell(number3);

				wbook.write();
				wbook.close();

			} catch (Exception e) {
				System.out.println(e);
			}
		}
	}
	

	public static void update_Infection(int firstPath, int path, int j)
	{
		if(firstPath!=path){
			for(int length=0;length<NODENUM;length++)
				if(PATH[firstPath].charAt(length)!=PATH[path].charAt(length)){
					if(PATH[firstPath].charAt(length)!=' '&&PATH[path].charAt(length)!=' ')
						infection[j][length]++;
				}
		}
		
	}

	// 根据status路径覆盖情况，轮盘赌选择一个一个未覆盖的路径，作为接下来的优化目标
	public static int random_UncoverPath(boolean[] status,int path,int repeat[]) {

		int[] similar = new int[PATHNUM];
		for(int i=0;i<PATHNUM;i++)
		{
			for(int j=0;j<NODENUM;j++)
			{
				if(PATH[path].charAt(j)==PATH[i].charAt(j))
					similar[i]++;			//统计path与所有路径之间的差异程度
			}
		}
		int[] possible = new int[PATHNUM];
		for(int i=0;i<PATHNUM;i++)
			if(!status[i])
				possible[i] = similar[i];
			else possible[i] = 0;
		//根据possible[]进行轮盘赌选择，目标路径
		int temp=0; int index=0;
		for(int i=0;i<PATHNUM;i++)
			temp+=possible[i]*Math.pow(0.5, repeat[i]);	//统计总数
		if (temp == 0) {
			index = (int)(Math.random()*PATHNUM);
			while(status[index]){
				index = (int)(Math.random()*PATHNUM);
			}
		}else{
			double rand[] = new double[PATHNUM];
			for(int i=0;i<PATHNUM;i++)
				rand[i]=((double)possible[i]*Math.pow(0.5, repeat[i]))/temp;		//统计各自的概率
			double random = Math.random();
			
			double bound[] = new double[PATHNUM];
			for(int i=0;i<PATHNUM;i++)
				for(int j=0;j<=i;j++)
					bound[i] += rand[j];	//计算轮盘赌刻度
			for(int i=0;i<PATHNUM;i++)
				if(random<bound[i]){
					index = i;
					break;
				}
		}
		return index;
	}
	
	public static int checksum(String ISBN){
		int sum=0; int k=0;
		for(int i=0;i<ISBN.length();i++){
			if(ISBN.charAt(i)-'0'>=0&&ISBN.charAt(i)-'0'<=9){
				k++;
				sum+=(ISBN.charAt(i)-'0')*k;
			}
			if(ISBN.charAt(i)=='X'||ISBN.charAt(i)=='x')
			{
				k++;
				sum+=10*k;
			}
		}
		return sum;
	}
	
	public static void DE_search(int[][]x, int[][]v, double[]fitness_x, double[]fitness_v, int[][] solution, boolean[] status, int run, int[] obj)
	{
		/* 差分变异操作 */
		int path;
       		for (int i = 0; i < pop_num; i++) {
			int k1 = (int) Math.floor(Math.random() * pop_num);
			while (k1 == i)
				k1 = (int) Math.floor(Math.random() * pop_num);
			int k2 = (int) Math.floor(Math.random() * pop_num);
			while (k2 == i || k2 == k1)
				k2 = (int) Math.floor(Math.random() * pop_num);
			int jrand = (int) (Math.random() * R);
			for (int j = 0; j < R; j++) {
				v[i][j] = (int) (x[i][j] + F * (x[k1][j] - x[k2][j]));
				if (Math.random() > Pc && j != jrand)
					v[i][j] = x[i][j];
				if (v[i][j] >= ub[j] || v[i][j] < lb[j]) {
					v[i][j] = (int) (lb[j] + Math.random() * (ub[j] - lb[j]));
				}
			}

			path = pathnum(v[i], fun_num);
			record[i] = path;
			if (!status[path]) {
				for (int j = 0; j < R; j++)
					solution[path][j] = v[i][j];
				status[path] = true;
				obj[0]++;
				nodeiscoverage(v[i], fun_num); // 标记已被覆盖的分支
			}

			if (obj[0] == PATHNUM)
				break;
		}

		for (int i = 0; i < pop_num; i++) {
			
			fitness_x[i] = benchmarkfunction(x[i], fun_num, -1);
			fitness_v[i] = benchmarkfunction(v[i], fun_num, -1);
			case_num[run] = case_num[run] + 2; // 评估次数更新

			if (fitness_v[i] > fitness_x[i]) // step 6：比较更新测试用例
			{
				for (int j = 0; j < R; j++)
					x[i][j] = v[i][j];
				fitness_x[i] = fitness_v[i];
			}
			if (obj[0] == PATHNUM) {
				break;
			}
		}
	}
	
	/*IPS版本的AVM算法*/
	public static void AVM_IPS(int[][]x, double[]fitness_x, int[][] solution, boolean[] status, int run, int[] obj)
	{
		//统计剩余路径被选为目标路径次数
		int repeat[] = new int [PATHNUM];
		for(int i=0;i<PATHNUM;i++)
			repeat[i] = 0;
		for(int i=0;i<pop_num;i++)  //遍历所有个体
		{
			if(obj[0] == PATHNUM)
				break;
			
			int target_path = random_UncoverPath(status,record[i],repeat);		//生成目标路径
			fitness_x[i] = benchmarkfunction(x[i], fun_num, target_path);
			int path;
			double[] fitness_temp = new double[2];
			repeat[target_path]++;
			
			int [][]temp = new int[2][R];
			for(int dim = 0;dim<R;dim++)
			{
				int step; //搜索步长
				boolean optimum = true;
				while(optimum)
				{
					//生成评估f(x-1)
					temp[0] = createTestCase(x[i], dim, x[i][dim]-1);
					if(temp[0][dim]>ub[dim]||temp[0][dim]<lb[dim])
					{
						fitness_temp[0] = -1;
					}else{
						path = pathnum(temp[0], fun_num);                       //获取覆盖的路径					
						fitness_temp[0] = benchmarkfunction(temp[0], fun_num, target_path);     //评估个体
						
						if(!status[path])
						{
							for(int t=0;t<R;t++)
								solution[path][t] = temp[0][t];                 //记录路径第一个
							status[path] = true;                             //标记路径Path是否已找到覆盖它的用例
							obj[0]++;                                           //已覆盖的路径数
							nodeiscoverage(temp[0],fun_num);                    //标记已被覆盖的分支
						}

						case_num[run] = case_num[run] + 1;           //评估次数更新
						if(case_num[run]>MCN)
							break;
						
						if(obj[0] == PATHNUM)
							break;
						
						if(status[target_path])					//如果选择的路径已经被覆盖，重新生成目标路径
							break;
					}
					
					
					
					//生成评估f(x+1)
					temp[1] = createTestCase(x[i], dim, x[i][dim]+1);
					if(temp[1][dim]>ub[dim]||temp[1][dim]<lb[dim])
					{
						fitness_temp[1] = -1;
					}else {
						path = pathnum(temp[1], fun_num);                       //获取覆盖的路径
						fitness_temp[1] = benchmarkfunction(temp[1], fun_num, target_path);     //评估个体
						
						if(!status[path])
						{
							for(int t=0;t<R;t++)
								solution[path][t] = temp[1][t];                 //记录路径第一个
							status[path] = true;                             //标记路径Path是否已找到覆盖它的用例
							obj[0]++;                                           //已覆盖的路径数
							nodeiscoverage(temp[1],fun_num);                    //标记已被覆盖的分支
						}
						
						case_num[run] = case_num[run] + 1;           //评估次数更新
						if(case_num[run]>MCN)
							break;
						
						if(obj[0] == PATHNUM)
							break;
						
						if(status[target_path])					//如果选择的路径已经被覆盖，重新生成目标路径
							break;
					}
					
					
					
					//第一层判断
					if(fitness_temp[0]<=fitness_x[i] && fitness_temp[1]<=fitness_x[i])
					{
						optimum = false;
						break;
					}
					
					//第二层判断
					if(fitness_temp[0]>fitness_temp[1]&&(temp[0][dim]<=ub[dim]&&temp[0][dim]>=lb[dim]))
						step = -1;
					else
						step = 1;
					
					//新增判断是否溢出
					if(fitness_temp[0]==-1&&fitness_temp[1]==-1)
					{
						optimum = false;
						break;
					}
					
					//第三层循环
					//生成评估f(x+step)
					while(benchmarkfunction(createTestCase(x[i], dim, x[i][dim]+step), fun_num, target_path)>fitness_x[i])
					{
						if(x[i][dim]>ub[dim]||x[i][dim]<lb[dim])
							break;
						x[i][dim] = x[i][dim] + step;
						
						path = pathnum(x[i], fun_num);                       //获取覆盖的路径
						fitness_x[i] = benchmarkfunction(x[i], fun_num, target_path);     //评估个体

						if(!status[path])
						{
							for(int t=0;t<R;t++)
								solution[path][t] = x[i][t];                 //记录路径第一个
							status[path] = true;                             //标记路径Path是否已找到覆盖它的用例
							obj[0]++;                                           //已覆盖的路径数
							nodeiscoverage(x[i],fun_num);                    //标记已被覆盖的分支
						}
						
						case_num[run] = case_num[run] + 2;           //评估次数更新（两次评估一次在循环判断条件，一次在更新x[i]）
						if(case_num[run]>MCN)
							break;
						
						if(obj[0] == PATHNUM)
							break;
						
						if(status[target_path])					//如果选择的路径已经被覆盖，重新生成目标路径
							break;		
						
						//更新搜索步长
						step = 2 * step;
						
						if(((x[i][dim]+step)>ub[dim]||(x[i][dim]+step)<lb[dim]))
							break;
					}
					if(case_num[run]>MCN)
						break;
				}
				if(case_num[run]>MCN)
					break;
			}
		}
	}
	
	/* 快速更新关系矩阵操作 */
	public static void Homo_Euclidean_Update(int[][]x, double[]fitness_x, int[][] solution, boolean[] status, int run, int[] obj)
	{
		//针对所有目标路径更新关系矩阵
		for(int p=0;p<PATHNUM;p++)
		{
			if(status[p]==true)
				break;
			int target_path = p;
			int path;
			int i = (int) (Math.random()*pop_num);
			//遍历所有维度，搜索出节点覆盖
			int dim = Math.random()>0.5?0:R-1;
			if(dim==0)
			{
				for(;dim<R;dim++)
				{
					int j = dim;                      //获取变更的变量下标
					int best;                                      //记录搜索出最优个体的值
					double[] fitness_temp = new double[step_length1+1];              //暂存变更后所有个体适应值
					
					int step;//初始化步长
					best = x[i][j];//初始化搜索最优个体值
					int firstPath = pathnum(x[i], fun_num);// 记录搜索前覆盖的路径
					
					step = (ub[j]-lb[j])/step_length1;
					
					while(step>1){//搜索个体数目大于step_length情况
						int[] temp = getIndex1(lb[j],ub[j],best,step);//temp暂存变更变量的值
						
						for(int k=0;k<step_length1+1;k++){
							x[i][j] = temp[k]; //替换搜索变量值
							if(x[i][j]>ub[j]||x[i][j]<lb[j])         //超出范围继续循环
								continue;
							
							path = pathnum(x[i], fun_num);                       //获取覆盖的路径
							update_Infection(firstPath, path, j);
							
							if(!status[path])
							{
								for(int t=0;t<R;t++)
									solution[path][t] = x[i][t];                 //记录路径第一个
								status[path] = true;                             //标记路径Path是否已找到覆盖它的用例
								obj[0]++;                                           //已覆盖的路径数
								nodeiscoverage(x[i],fun_num);                    //标记已被覆盖的分支
							}
							
							if(obj[0] == PATHNUM)
								break;
							
							if(status[target_path])					//如果选择的路径已经被覆盖，重新生成目标路径
								break;
							
							fitness_temp[k] = benchmarkfunction(x[i], fun_num, target_path);     //评估个体
							case_num[run] = case_num[run] + 1;           //评估次数更新
						}
						int best_index = getBestIndex(fitness_temp);
						x[i][j] = temp[best_index];
						best = temp[best_index];
						fitness_x[i] = fitness_temp[best_index];
						
						step = step/step_length1;                                  //step下降更新
					}
					//搜索个体数目小于step_length情况
					step = 1;
					int[] temp = getIndex1(lb[j],ub[j],best,step);
					
					for(int k=0;k<step_length1+1;k++)
					{
						x[i][j] = temp[k]; // 替换搜索变量值
						if(x[i][j]>ub[j]||x[i][j]<lb[j])         //超出范围继续循环
							continue;
						
						path = pathnum(x[i], fun_num); // 获取覆盖的路径
						update_Infection(firstPath, path, j);
						
						if (!status[path]) {
							for (int t = 0; t < R; t++)
								solution[path][t] = x[i][t]; // 记录路径第一个
							status[path] = true; // 标记路径Path是否已找到覆盖它的用例
							obj[0]++; // 已覆盖的路径数
							nodeiscoverage(x[i], fun_num); // 标记已被覆盖的分支
						}
						if (obj[0] == PATHNUM)
							break;
						
						if(status[target_path])					//如果选择的路径已经被覆盖，重新生成目标路径
							break;

						fitness_temp[k] = benchmarkfunction(x[i], fun_num, target_path); // 评估个体
						case_num[run] = case_num[run] + 1; // 评估次数更新
					}
					
					int best_index = getBestIndex(fitness_temp);
					x[i][j] = temp[best_index];
					fitness_x[i] = fitness_temp[best_index];
				}
			}else{
				for(;dim>=0;dim--)
				{
					int j = dim;                      //获取变更的变量下标
					int best;                                      //记录搜索出最优个体的值
					double[] fitness_temp = new double[step_length1+1];              //暂存变更后所有个体适应值
					
					int step;//初始化步长
					best = x[i][j];//初始化搜索最优个体值
					int firstPath = pathnum(x[i], fun_num);// 记录搜索前覆盖的路径
					
					step = (ub[j]-lb[j])/step_length1;
					
					while(step>1){//搜索个体数目大于step_length情况
						int[] temp = getIndex1(lb[j],ub[j],best,step);//temp暂存变更变量的值
						
						for(int k=0;k<step_length1+1;k++){
							x[i][j] = temp[k]; //替换搜索变量值
							if(x[i][j]>ub[j]||x[i][j]<lb[j])         //超出范围继续循环
								continue;
							
							path = pathnum(x[i], fun_num);                       //获取覆盖的路径
							update_Infection(firstPath, path, j);
							
							if(!status[path])
							{
								for(int t=0;t<R;t++)
									solution[path][t] = x[i][t];                 //记录路径第一个
								status[path] = true;                             //标记路径Path是否已找到覆盖它的用例
								obj[0]++;                                           //已覆盖的路径数
								nodeiscoverage(x[i],fun_num);                    //标记已被覆盖的分支
							}
							
							if(obj[0] == PATHNUM)
								break;
							
							if(status[target_path])					//如果选择的路径已经被覆盖，重新生成目标路径
								break;
							
							fitness_temp[k] = benchmarkfunction(x[i], fun_num, target_path);     //评估个体
							case_num[run] = case_num[run] + 1;           //评估次数更新
						}
						int best_index = getBestIndex(fitness_temp);
						x[i][j] = temp[best_index];
						best = temp[best_index];
						fitness_x[i] = fitness_temp[best_index];
						
						step = step/step_length1;                                  //step下降更新
					}
					//搜索个体数目小于step_length情况
					step = 1;
					int[] temp = getIndex1(lb[j],ub[j],best,step);
					
					for(int k=0;k<step_length1+1;k++)
					{
						x[i][j] = temp[k]; // 替换搜索变量值
						if(x[i][j]>ub[j]||x[i][j]<lb[j])         //超出范围继续循环
							continue;
						
						path = pathnum(x[i], fun_num); // 获取覆盖的路径
						update_Infection(firstPath, path, j);
						
						if (!status[path]) {
							for (int t = 0; t < R; t++)
								solution[path][t] = x[i][t]; // 记录路径第一个
							status[path] = true; // 标记路径Path是否已找到覆盖它的用例
							obj[0]++; // 已覆盖的路径数
							nodeiscoverage(x[i], fun_num); // 标记已被覆盖的分支
						}
						if (obj[0] == PATHNUM)
							break;
						
						if(status[target_path])					//如果选择的路径已经被覆盖，重新生成目标路径
							break;

						fitness_temp[k] = benchmarkfunction(x[i], fun_num, target_path); // 评估个体
						case_num[run] = case_num[run] + 1; // 评估次数更新
					}
					
					int best_index = getBestIndex(fitness_temp);
					x[i][j] = temp[best_index];
					fitness_x[i] = fitness_temp[best_index];
				}
			}
			
			if(obj[0] == PATHNUM)
				break;                                   //判断路径是否全部覆盖，如果全部覆盖则退出循环
		}
			
	}
	
	public static void ScatterSearch(int[][]x, double[]fitness_x, int[][] solution, boolean[] status, int run, int[] obj)
	{
		int repeat[] = new int [PATHNUM];
		for(int i=0;i<PATHNUM;i++)
			repeat[i] = 0;
		for(int i=0;i<pop_num;i++)					//遍历所有个体
		{
			/* 根据个体覆盖路径与剩余路径相似程度（由走向相同节点数目评估），轮盘赌的形式选择一条路径作为目标。
			 * 接下来RP搜索过程就是使得优化个体，向着选择的目标路径方法搜索*/
			int target_path = random_UncoverPath(status,record[i],repeat);		//生成目标路径
			int path;
			repeat[target_path]++;
			fitness_x[i] = benchmarkfunction(x[i], fun_num, target_path);
			
			int[][]temp_x = new int[step_length1][R];
			
			//遍历所有维度，搜索出节点覆盖
			int dim = Math.random()>0.5?0:R-1;
			if(dim==0)
			{
				for(;dim<R;dim++)
				{
					int j = dim;                      //获取变更的变量下标
					int best;                                      //记录搜索出最优个体的值
					double[] fitness_temp = new double[step_length1];              //暂存变更后所有个体适应值
					
					int step;//初始化步长
					best = x[i][j];//初始化搜索最优个体值
					int firstPath = pathnum(x[i], fun_num);// 记录搜索前覆盖的路径
					
					step = (ub[j]-lb[j])/step_length1;
					
					while(step>1)
					{//搜索个体数目大于step_length情况
						int[] temp = getIndex3(lb[j],ub[j],best,step);//temp暂存变更变量的值
						
						for(int k=0;k<step_length1;k++)
						{
//							x[i][j] = temp[k]; //替换搜索变量值
							temp_x[k] = createTestCase(x[i], j, temp[k]);
							if(temp_x[k][j]>ub[j]||temp_x[k][j]<lb[j])         //超出范围继续循环
								continue;
							
							path = pathnum(temp_x[k], fun_num);                       //获取覆盖的路径
							update_Infection(firstPath, path, j);
							fitness_temp[k] = benchmarkfunction(temp_x[k], fun_num, target_path);     //评估个体
							case_num[run] = case_num[run] + 1;           //评估次数更新
							
							if(!status[path])
							{
								for(int t=0;t<R;t++)
									solution[path][t] = temp_x[k][t];                 //记录路径第一个
								status[path] = true;                             //标记路径Path是否已找到覆盖它的用例
								obj[0]++;                                           //已覆盖的路径数
								nodeiscoverage(temp_x[k],fun_num);                    //标记已被覆盖的分支
							}
							
							if(obj[0] == PATHNUM)
								break;
							
							if(status[target_path])					//如果选择的路径已经被覆盖，重新生成目标路径
								break;
						}
						int best_index = getBestIndex(fitness_temp);
						if(fitness_temp[best_index]>fitness_x[i])
						{
							x[i][j] = temp[best_index];
							best = temp[best_index];
							fitness_x[i] = fitness_temp[best_index];
						}	
						step = step/step_length1;                                  //step下降更新
					}
					//搜索个体数目小于step_length情况
					step = 1;
					int[] temp = getIndex3(lb[j],ub[j],best,step);
					
					for(int k=0;k<step_length1;k++)
					{
//						x[i][j] = temp[k]; // 替换搜索变量值
						temp_x[k] = createTestCase(x[i], j, temp[k]);
						if(temp_x[k][j]>ub[j]||temp_x[k][j]<lb[j])         //超出范围继续循环
							continue;
						
						path = pathnum(temp_x[k], fun_num);                       //获取覆盖的路径
						update_Infection(firstPath, path, j);
						fitness_temp[k] = benchmarkfunction(temp_x[k], fun_num, target_path);     //评估个体
						case_num[run] = case_num[run] + 1;           //评估次数更新
						
						if(!status[path])
						{
							for(int t=0;t<R;t++)
								solution[path][t] = temp_x[k][t];                 //记录路径第一个
							status[path] = true;                             //标记路径Path是否已找到覆盖它的用例
							obj[0]++;                                           //已覆盖的路径数
							nodeiscoverage(temp_x[k],fun_num);                    //标记已被覆盖的分支
						}
						
						if(obj[0] == PATHNUM)
							break;
						
						if(status[target_path])					//如果选择的路径已经被覆盖，重新生成目标路径
							break;
					}
					int best_index = getBestIndex(fitness_temp);
					if(fitness_temp[best_index]>fitness_x[i])
					{
						x[i][j] = temp[best_index];
						best = temp[best_index];
						fitness_x[i] = fitness_temp[best_index];
					}
				}
			}else{
				for(;dim<R;dim++)
				{
					int j = dim;                      //获取变更的变量下标
					int best;                                      //记录搜索出最优个体的值
					double[] fitness_temp = new double[step_length1];              //暂存变更后所有个体适应值
					
					int step;//初始化步长
					best = x[i][j];//初始化搜索最优个体值
					int firstPath = pathnum(x[i], fun_num);// 记录搜索前覆盖的路径
					
					step = (ub[j]-lb[j])/step_length1;
					
					while(step>1)
					{//搜索个体数目大于step_length情况
						int[] temp = getIndex3(lb[j],ub[j],best,step);//temp暂存变更变量的值
						
						for(int k=0;k<step_length1;k++)
						{
//							x[i][j] = temp[k]; //替换搜索变量值
							temp_x[k] = createTestCase(x[i], j, temp[k]);
							if(temp_x[k][j]>ub[j]||temp_x[k][j]<lb[j])         //超出范围继续循环
								continue;
							
							path = pathnum(temp_x[k], fun_num);                       //获取覆盖的路径
							update_Infection(firstPath, path, j);
							fitness_temp[k] = benchmarkfunction(temp_x[k], fun_num, target_path);     //评估个体
							case_num[run] = case_num[run] + 1;           //评估次数更新
							
							if(!status[path])
							{
								for(int t=0;t<R;t++)
									solution[path][t] = temp_x[k][t];                 //记录路径第一个
								status[path] = true;                             //标记路径Path是否已找到覆盖它的用例
								obj[0]++;                                           //已覆盖的路径数
								nodeiscoverage(temp_x[k],fun_num);                    //标记已被覆盖的分支
							}
							
							if(obj[0] == PATHNUM)
								break;
							
							if(status[target_path])					//如果选择的路径已经被覆盖，重新生成目标路径
								break;
						}
						int best_index = getBestIndex(fitness_temp);
						if(fitness_temp[best_index]>fitness_x[i])
						{
							x[i][j] = temp[best_index];
							best = temp[best_index];
							fitness_x[i] = fitness_temp[best_index];
						}	
						step = step/step_length1;                                  //step下降更新
					}
					//搜索个体数目小于step_length情况
					step = 1;
					int[] temp = getIndex3(lb[j],ub[j],best,step);
					
					for(int k=0;k<step_length1;k++)
					{
//						x[i][j] = temp[k]; // 替换搜索变量值
						temp_x[k] = createTestCase(x[i], j, temp[k]);
						if(temp_x[k][j]>ub[j]||temp_x[k][j]<lb[j])         //超出范围继续循环
							continue;
						
						path = pathnum(temp_x[k], fun_num);                       //获取覆盖的路径
						update_Infection(firstPath, path, j);
						fitness_temp[k] = benchmarkfunction(temp_x[k], fun_num, target_path);     //评估个体
						case_num[run] = case_num[run] + 1;           //评估次数更新
						
						if(!status[path])
						{
							for(int t=0;t<R;t++)
								solution[path][t] = temp_x[k][t];                 //记录路径第一个
							status[path] = true;                             //标记路径Path是否已找到覆盖它的用例
							obj[0]++;                                           //已覆盖的路径数
							nodeiscoverage(temp_x[k],fun_num);                    //标记已被覆盖的分支
						}
						
						if(obj[0] == PATHNUM)
							break;
						
						if(status[target_path])					//如果选择的路径已经被覆盖，重新生成目标路径
							break;
					}
					int best_index = getBestIndex(fitness_temp);
					if(fitness_temp[best_index]>fitness_x[i])
					{
						x[i][j] = temp[best_index];
						best = temp[best_index];
						fitness_x[i] = fitness_temp[best_index];
					}
				}
			}
			if(obj[0] == PATHNUM)
				break;                                   //判断路径是否全部覆盖，如果全部覆盖则退出循环
		}
	}
	
	public static void RelationshipProgramming(int[][]x, double[]fitness_x, int[][] solution, boolean[] status, int run, int[] obj)
	{
		int repeat[] = new int[PATHNUM];
		for(int i=0;i<PATHNUM;i++)
			repeat[i] = 0;
		for (int i = 0; i < pop_num; i++) // 遍历所有个体
		{
			/*
			 * 根据个体覆盖路径与剩余路径相似程度（由走向相同节点数目评估），轮盘赌的形式选择一条路径作为目标。
			 * 接下来RP搜索过程就是使得优化个体，向着选择的目标路径方法搜索
			 */
			int target_path = random_UncoverPath(status, record[i],repeat); // 生成目标路径
			int path;
			repeat[target_path]++;
			
			for (int node = 0; node < PATH[target_path].length(); node++) {
				if (PATH[target_path].charAt(node) != PATH[record[i]].charAt(node)) {
					// 如果需要覆盖的路径和测试用例当前覆盖的路径的分支节点走向对比不一样，则对影响节点走向维度进行搜索
					int tempare = 0;
					for (int r = 0; r < R; r++)
						tempare += infection[r][node];

					int firstPath; // 记录初始测试用例覆盖的路径

					if (tempare == 0) {// 如果没有维度影响，此时只能随机
						int j = (int) (Math.random() * R); // 随机获取变更的变量下标
						
						int best; // 记录搜索出最优个体的值
						double[] fitness_temp = new double[step_length2 + 1]; // 暂存变更后所有个体适应值

						int step;// 初始化步长
						best = x[i][j];// 初始化搜索最优个体值
						firstPath = pathnum(x[i], fun_num);// 记录搜索前覆盖的路径

						step = (ub[j] - lb[j]) / step_length2;

						while (step > 1) {// 搜索个体数目大于step_length情况
							int[] temp = getIndex2(lb[j], ub[j], best, step);// temp暂存变更变量的值

							for (int k = 0; k < step_length2 + 1; k++) {
								x[i][j] = temp[k]; // 替换搜索变量值
								if(x[i][j]>ub[j]||x[i][j]<lb[j])         //超出范围继续循环
									continue;
								
								path = pathnum(x[i], fun_num); // 获取覆盖的路径

								update_Infection(firstPath, path, j);

								if (!status[path]) {
									for (int t = 0; t < R; t++)
										solution[path][t] = x[i][t]; // 记录路径第一个
									status[path] = true; // 标记路径Path是否已找到覆盖它的用例
									obj[0]++; // 已覆盖的路径数
								}

								if (obj[0] == PATHNUM)
									break;

								if (status[target_path]) // 如果选择的路径已经被覆盖，重新生成目标路径
									break;

								fitness_temp[k] = benchmarkfunction(x[i], fun_num, target_path); // 评估个体
								case_num[run] = case_num[run] + 1; // 评估次数更新
								
//								System.out.println(1);
							}
							int best_index = getBestIndex(fitness_temp);
							x[i][j] = temp[best_index];
							best = temp[best_index];
							fitness_x[i] = fitness_temp[best_index];

							step = step / step_length2; // step下降更新

							if (PATH[target_path].charAt(node) == PATH[record[i]].charAt(node))
								break;
						}
						// 搜索个体数目小于step_length情况
						step = 1;
						int[] temp = getIndex2(lb[j], ub[j], best, step);
						for (int k = 0; k < step_length2 + 1; k++) {
							x[i][j] = temp[k]; // 替换搜索变量值
							if(x[i][j]>ub[j]||x[i][j]<lb[j])         //超出范围继续循环
								continue;
							
							path = pathnum(x[i], fun_num); // 获取覆盖的路径

							/*
							 * 判断在维度搜索过程中算法是否有覆盖到一条不同的路径，
							 * 如果有，则说明该维度的改变对节点的走向有影响
							 */
							update_Infection(firstPath, path, j);

							if (!status[path]) {
								for (int t = 0; t < R; t++)
									solution[path][t] = x[i][t]; // 记录路径第一个
								status[path] = true; // 标记路径Path是否已找到覆盖它的用例
								obj[0]++; // 已覆盖的路径数
							}
							if (obj[0] == PATHNUM)
								break;

							if (status[target_path]) // 如果选择的路径已经被覆盖，重新生成目标路径
								break;

							fitness_temp[k] = benchmarkfunction(x[i], fun_num, target_path); // 评估个体
							case_num[run] = case_num[run] + 1; // 评估次数更新
							
//							System.out.println(2);
						}

						int best_index = getBestIndex(fitness_temp);
						x[i][j] = temp[best_index];
						fitness_x[i] = fitness_temp[best_index];

					} else {// 有部分或者全部维度有影响，则进行轮盘赌选择维度优化
						double rand[] = new double[R];
						for (int r = 0; r < R; r++)
							rand[r] = ((double) infection[r][node] + 1) / (tempare + R);
						// 由相关度，轮盘赌选择维度
						double random = Math.random();
						int index = 0;
						double bound[] = new double[R];
						for (int r = 0; r < R; r++)
							for (int j = 0; j <= r; j++)
								bound[r] += rand[j]; // 计算轮盘赌刻度
						for (int r = 0; r < R; r++)
							if (random < bound[r]) {
								index = r; // 获取轮盘赌赌博得到的刻度，选择维度index进行优化
								break;
							}
						int j = index;

						firstPath = pathnum(x[i], fun_num);// 记录搜索前覆盖的路径

						/* 单维度搜索：先大步长，后小步长 */
						int best; // 记录搜索出最优个体的值
						double[] fitness_temp = new double[step_length2 + 1]; // 暂存变更后所有个体适应值

						double step;// 初始化步长
						best = x[i][index];// 初始化搜索最优个体值
						step = (ub[index] - lb[index]) / step_length2;

						while (step > 1) {// 搜索个体数目大于step_length情况
							int[] temp = getIndex2(lb[j], ub[j], best, step);// temp暂存变更变量的值

							for (int k = 0; k < step_length2 + 1; k++) {
								x[i][j] = temp[k]; // 替换搜索变量值
								if(x[i][j]>ub[j]||x[i][j]<lb[j])         //超出范围继续循环
									continue;
								
								path = pathnum(x[i], fun_num); // 获取覆盖的路径

								update_Infection(firstPath, path, j);

								if (!status[path]) {
									for (int t = 0; t < R; t++)
										solution[path][t] = x[i][t]; // 记录路径第一个
									status[path] = true; // 标记路径Path是否已找到覆盖它的用例
									obj[0]++; // 已覆盖的路径数
								}

								if (obj[0] == PATHNUM)
									break;

								if (status[target_path]) // 如果选择的路径已经被覆盖，重新生成目标路径
									break;

								fitness_temp[k] = benchmarkfunction(x[i], fun_num, target_path); // 评估个体
								case_num[run] = case_num[run] + 1; // 评估次数更新
								
//								System.out.println(3);
							}
							int best_index = getBestIndex(fitness_temp);
							x[i][j] = temp[best_index];
							best = temp[best_index];
							fitness_x[i] = fitness_temp[best_index];

							step = step / step_length2; // step下降更新

							if (PATH[target_path].charAt(node) == PATH[record[i]].charAt(node))
								break;
						}
						// 搜索个体数目小于step_length情况
						step = 1;
						int[] temp = getIndex2(lb[j], ub[j], best, step);
						for (int k = 0; k < step_length2 + 1; k++) 
						{
							x[i][j] = temp[k]; // 替换搜索变量值
							if(x[i][j]>ub[j]||x[i][j]<lb[j])         //超出范围继续循环
								continue;
							
							path = pathnum(x[i], fun_num); // 获取覆盖的路径

							/*
							 * 判断在维度搜索过程中算法是否有覆盖到一条不同的路径，
							 * 如果有，则说明该维度的改变对节点的走向有影响
							 */
							update_Infection(firstPath, path, j);

							if (!status[path]) {
								for (int t = 0; t < R; t++)
									solution[path][t] = x[i][t]; // 记录路径第一个
								status[path] = true; // 标记路径Path是否已找到覆盖它的用例
								obj[0]++; // 已覆盖的路径数
							}
							if (obj[0] == PATHNUM)
								break;

							if (status[target_path]) // 如果选择的路径已经被覆盖，重新生成目标路径
								break;

							fitness_temp[k] = benchmarkfunction(x[i], fun_num, target_path); // 评估个体
							case_num[run] = case_num[run] + 1; // 评估次数更新
							
//							System.out.println(case_num[run]);
						}

						int best_index = getBestIndex(fitness_temp);
						x[i][j] = temp[best_index];
						fitness_x[i] = fitness_temp[best_index];
					}
				}
			}
			if (obj[0] == PATHNUM)
				break; // 判断路径是否全部覆盖，如果全部覆盖则退出循环
		}
	}
	
	/*
	 * 3个benchmark函数：1)transmit ;2)send ;3)processTupleArrival ;4)executeTuple ;5)checkCloudletCompletion ;6)getRsultantTuples
	 * (函数序号，输入变量维度，存在路径数目，可行路径数目)
	 * transmit(1,3,2,1) 100%; send(2,2,9,5)66%; processEvent(3,7,9,6)100%; executeTuple(4,7,5,3)100%; 
	 * checkCloudletCompletion(5,5,6,3)100%; getResultantTuple(6,8,7,4)100%
	 * getPhysicalTopology(7,110,7,6);  
	 * gimp_rgb_to_hsv(7,3,6,5);  gimp_hsv_to_rgb(8,3,15,3)
	 * check_ISBN(9,12,7,5);  check_ISSN(10,10,7,5);  gimp_hwb_to_rgb(11,3,15,3)
	 * TIFF_GetSourceSamples(12,2,8,7)
	 */
	public static void setFunctionParameters(){
		if (fun_num == 1) {
			R = 3;
			PATHNUM = 2;
			NODENUM = 1;
			col = 1;
		}
		if (fun_num == 2) {
			R = 2;
			PATHNUM = 9;
			NODENUM = 5;
			col = 2;
		}
		if (fun_num == 3) {
			R = 7;
			PATHNUM = 9;
			NODENUM = 6;
			col = 3;
		}
		if (fun_num == 4) {
			R = 7;
			PATHNUM = 5;
			NODENUM = 3;
			col = 4;
		}
		if (fun_num == 5) {
			R = 5;
			PATHNUM = 6;
			NODENUM = 3;
			col = 5;
		}
		if (fun_num == 6) {
			R = 8;
			PATHNUM = 7;
			NODENUM = 4;
			col = 6;
		}
//		if (fun_num == 7) {
//			R = 110;
//			PATHNUM = 7;
//			NODENUM = 6;
//			col = 7;
//		}
		if (fun_num == 7) {
			R = 3;
			PATHNUM = 6;
			NODENUM = 5;
			col = 7;
		}
		if (fun_num == 8) {
			R = 3;
			PATHNUM = 15;
			NODENUM = 3;
			col = 8;
		}
		if (fun_num == 9) {
			R = 12;
			PATHNUM = 7;
			NODENUM = 5;
			col = 9;
		}
		if (fun_num == 10) {
			R = 10;
			PATHNUM = 7;
			NODENUM = 5;
			col = 10;
		}
		if (fun_num == 11) {
			R = 3;
			PATHNUM = 15;
			NODENUM = 3;
			col = 11;
		}
		if (fun_num == 12) {
			R = 2;
			PATHNUM = 8;
			NODENUM = 7;
			col = 12;
		}
		if (fun_num == 13) {
			R = 7;
			PATHNUM = 48;
			NODENUM = 4;
			col = 13;
		}
		if (fun_num == 14) {
			R = 6;
			PATHNUM = 3;
			NODENUM = 2;
			col = 14;
		}
		if (fun_num == 15) {
			R = 11;
			PATHNUM = 12;
			NODENUM = 7;
			col = 15;
		}
		if (fun_num == 16) {
			R = 4;
			PATHNUM = 3;
			NODENUM = 2;
			col = 16;
		}
		if (fun_num == 17) {
			R = 8;
			PATHNUM = 4;
			NODENUM = 2;
			col = 17;
		}
		if (fun_num == 18) {
			R = 6;
			PATHNUM = 10;
			NODENUM = 5;
			col = 18;
		}
		if (fun_num == 19) {
			R = 3;
			PATHNUM = 4;
			NODENUM = 4;
			col = 19;
		}
		if (fun_num == 20) {
			R = 1;
			PATHNUM = 2;
			NODENUM = 1;
			col = 20;
		}
		if (fun_num == 21) {
			R = 10;
			PATHNUM = 2;
			NODENUM = 1;
			col = 21;
		}
		if (fun_num == 22) {
			R = 2;
			PATHNUM = 4;
			NODENUM = 2;
			col = 22;
		}
		if (fun_num == 23) {
			R = 3;
			PATHNUM = 4;
			NODENUM = 3;
			col = 23;
		}
		if (fun_num == 24) {
			R = 4;
			PATHNUM = 10;
			NODENUM = 5;
			col = 24;
		}
		if (fun_num == 25) {
			R = 3;
			PATHNUM = 15;
			NODENUM = 1;
			col = 25;
		}
		if (fun_num == 26) {
			R = 3;
			PATHNUM = 3;
			NODENUM = 2;
			col = 26;
		}
		if (fun_num == 27) {
			R = 2;
			PATHNUM = 11;
			NODENUM = 10;
			col = 27;
		}
		if (fun_num == 28) {
			R = 4;
			PATHNUM = 11;
			NODENUM = 4;
			col = 28;
		}
	}
	
	public static void setDomainAndEncoding() {
		// transmit(1,3,2,1) 100%;
		if (fun_num == 1) {
			for (int i = 0; i < R; i++) {
				lb[i] = 0; // 初始化上下界
				ub[i] = 255;
			}
		} // send(2,2,9,5)66%;
		if (fun_num == 2) {
			for (int i = 0; i < R; i++) {
				lb[i] = -1000000000; // 初始化上下界
				ub[i] = 1000000000;
			}
		} // processEvent(3,7,9,6)100%;
		if (fun_num == 3) {
			lb[0] = 0;
			ub[0] = Integer.MAX_VALUE;
			lb[1] = 0;
			ub[1] = 4;
			lb[2] = -2;
			ub[2] = Integer.MAX_VALUE-2;
			lb[3] = 0;
			ub[3] = Integer.MAX_VALUE;
			lb[4] = 0;
			ub[4] = Integer.MAX_VALUE;
			lb[5] = 0;
			ub[5] = Integer.MAX_VALUE;
			lb[6] = -2;
			ub[6] = Integer.MAX_VALUE-2;
		} // executeTuple(4,7,5,3)100%;
		if (fun_num == 4) {
			lb[0] = 1;
			ub[0] = 3;
			for (int i = 1; i < R; i++) {
				lb[i] = 0; // 初始化上下界
				ub[i] = 255;
			}
		} // checkCloudletCompletion(5,5,6,3)100%;
		if (fun_num == 5) {
			lb[0] = 0;
			ub[0] = 2;
			for (int i = 1; i < R - 1; i++) {
				lb[i] = 0; // 初始化上下界
				ub[i] = 255;
			}
			lb[R - 1] = 0;
			ub[R - 1] = 2;
		} // getResultantTuple(6,8,7,4)100%*/
		if (fun_num == 6) {
			for (int i = 0; i < 6; i++) {
				lb[i] = 0; // 初始化上下界
				ub[i] = 255;
			}
			lb[6] = 0;
			ub[6] = 2;
			lb[7] = 1;
			lb[7] = 4;
		} // getPhysicalTopology(7,110,7,6)
//		if (fun_num == 7) {
//			for (int i = 0; i < 10; i++) {
//				for (int j = 0; j < 10; j++) {
//					lb[j + i * 11] = 0;
//					ub[j + i * 11] = 255;
//				}
//				lb[10 + i * 11] = 0;
//				ub[10 + i * 11] = 4;
//			}
//		} // gimp_rgb_to_hsv(7,3,6,5);
		if (fun_num == 7) {
			for (int i = 0; i < R; i++) {
				lb[i] = 0;
				ub[i] = Integer.MAX_VALUE;
			}
		} // gimp_hsv_to_rgb(8,3,15,3)
		if (fun_num == 8) {
			lb[0] = 0;
			ub[0] = 1;
			for (int i = 1; i < R; i++) {
				lb[i] = 0;
				ub[i] = Integer.MAX_VALUE;
			}
		} // check_ISBN(9,12,7,5)
		if (fun_num == 9) {
			for (int i = 0; i < R - 1; i++) {
				lb[i] = 0;
//				ub[i] = 255;
				ub[i] = Integer.MAX_VALUE;
			}
			lb[R - 1] = 0;
			ub[R - 1] = 12;
		} // check_ISSN(10,10,7,5)
		if (fun_num == 10) {
			for (int i = 0; i < R - 1; i++) {
				lb[i] = 0;
				ub[i] = Integer.MAX_VALUE;
			}
			lb[R - 1] = 0;
			ub[R - 1] = 10;
		} // gimp_hwb_to_rgb(11,3,15,3)
		if (fun_num == 11) {
			lb[0] = -60;
			ub[0] = 365;
			lb[1] = 0;
			ub[1] = 1;
			lb[2] = 0;
			ub[2] = 1;
		} // TIFF_GetSourceSamples(12,2,8,7)
		if (fun_num == 12) {
			lb[0] = 0;
			ub[0] = Integer.MAX_VALUE;
			lb[1] = 0;
			ub[1] = Integer.MAX_VALUE;
		}
		// initFactory(13,7,48,4)
		if (fun_num == 13) {
			for (int i = 0; i < 6; i++) {
				lb[i] = 0;
				ub[i] = Integer.MAX_VALUE;
			}
			lb[6] = 0;
			ub[6] = Integer.MAX_VALUE;
			// ub[6] = 100;
		} // CleanXmlAnnotator(14,6,3,2)
		if (fun_num == 14) {
			for (int i = 0; i < R; i++) {
				lb[i] = 0; // 初始化上下界
				ub[i] = Integer.MAX_VALUE;
			}
		} // WordsToSentencesAnnotator(15,11,12,7)
		if (fun_num == 15) {
			lb[0] = 0;
			ub[0] = 1;
			lb[1] = 0;
			ub[1] = 1;
			for (int i = 2; i < R; i++) {
				lb[i] = 0; // 初始化上下界
				ub[i] = Integer.MAX_VALUE;
			}
		} // annotate(16,4,3,2)
		if (fun_num == 16) {
			for (int i = 0; i < R - 1; i++) {
				lb[i] = 0; // 初始化上下界
				ub[i] = Integer.MAX_VALUE;
			}
			lb[R - 1] = 1;
			ub[R - 1] = Integer.MAX_VALUE;
		} // NERClassifierCombiner(17,8,4,2)
		if (fun_num == 17) {
			for (int i = 0; i < R - 1; i++) {
				lb[i] = 0; // 初始化上下界
				ub[i] = Integer.MAX_VALUE;
			}
			lb[R - 1] = 0;
			ub[R - 1] = 1;
		} // setTrueCaseText(18,6,10,5)
		if (fun_num == 18) {
			for (int i = 0; i < R - 1; i++) {
				lb[i] = 0; // 初始化上下界
				ub[i] = Integer.MAX_VALUE;
			}
			lb[R - 1] = 0;
			ub[R - 1] = 1;
		} // Triangle(3,4,4);Factorial(1,2,1);BubbleSort(10,2,1);GCD(2,4,2);Middle(3,4,3);Commission(3,3,2);decision(4,11,4)
		if( (fun_num == 19) || (fun_num == 20) || (fun_num == 21) || (fun_num == 22) || (fun_num == 23) || (fun_num == 26) || (fun_num == 28))
			for(int j = 0 ; j < R ; j++ )
		    {
			    lb[j] = 1 ;
			    ub[j] = Integer.MAX_VALUE ;
		    }
		// Tomorrow(2,6,3)
		if(fun_num == 24)
		{
			lb[0] = 1 ;
			ub[0] = 7 ;
			lb[1] = 1900 ;
			ub[1] = Integer.MAX_VALUE ;
			lb[2] = 1 ;
			ub[2] = 12 ;
			lb[3] = 1 ;
			ub[3] = 31;
		} // calculator(2,15,15)
		if(fun_num == 25)
			for(int j = 0 ; j < R ; j++ )
		    {
			    lb[j] = 1 ;
			    ub[j] = 256 ;
		    }
		// Premium(2,11,10)
		if(fun_num == 27)
		{
			lb[0] = 1 ;
			ub[0] = Integer.MAX_VALUE ;
			lb[1] = 1 ;
			ub[1] = 12 ;
		}

		switch (fun_num) {
		case 1:
			PATH[0] = "0";
			PATH[1] = "1";
			break;
		case 2:
			PATH[0] = "0    ";
			PATH[1] = "100  ";
			PATH[2] = "1010 ";
			PATH[3] = "10110";
			PATH[4] = "10111";
			PATH[5] = "110  ";
			PATH[6] = "1110 ";
			PATH[7] = "11110";
			PATH[8] = "11111";
			break;
		case 3:
			PATH[0] = "0     ";
			PATH[1] = "10    ";
			PATH[2] = "110   ";
			PATH[3] = "111 00";
			PATH[4] = "111 01";
			PATH[5] = "111 1 ";
			PATH[6] = "12 0  ";
			PATH[7] = "12 1  ";
			PATH[8] = "13    ";
			break;
		case 4:
			PATH[0] = "000";
			PATH[1] = "001";
			PATH[2] = "010";
			PATH[3] = "011";
			PATH[4] = "1  ";
			break;
		case 5:
			PATH[0] = "000";
			PATH[1] = "001";
			PATH[2] = "010";
			PATH[3] = "011";
			PATH[4] = "1 0";
			PATH[5] = "1 1";
			break;
		case 6:
			PATH[0] = "0000";
			PATH[1] = "0001";
			PATH[2] = "001 ";
			PATH[3] = "0100";
			PATH[4] = "0101";
			PATH[5] = "011 ";
			PATH[6] = "1   ";
			break;
//		case 7:
//			PATH[0] = "0     ";
//			PATH[1] = "100   ";
//			PATH[2] = "1010  ";
//			PATH[3] = "10110 ";
//			PATH[4] = "10111 ";
//			PATH[5] = "11   0";
//			PATH[6] = "11   1";
//			break;
		case 7:
			PATH[0] = "000  ";
			PATH[1] = "001  ";
			PATH[2] = "01 0 ";
			PATH[3] = "01 10";
			PATH[4] = "01 11";
			PATH[5] = "1    ";
			break;
		case 8:
			PATH[0] = "0  ";
			PATH[1] = "100";
			PATH[2] = "101";
			PATH[3] = "102";
			PATH[4] = "103";
			PATH[5] = "104";
			PATH[6] = "105";
			PATH[7] = "106";
			PATH[8] = "110";
			PATH[9] = "111";
			PATH[10] = "112";
			PATH[11] = "113";
			PATH[12] = "114";
			PATH[13] = "115";
			PATH[14] = "116";
			break;
		case 9:
		case 10:
			PATH[0] = "0    ";
			PATH[1] = "10   ";
			PATH[2] = "1100 ";
			PATH[3] = "1101 ";
			PATH[4] = "111  ";
			PATH[5] = "2   0";
			PATH[6] = "2   1";
			break;
		case 11:
			PATH[0] = "0  ";
			PATH[1] = "100";
			PATH[2] = "101";
			PATH[3] = "102";
			PATH[4] = "103";
			PATH[5] = "104";
			PATH[6] = "105";
			PATH[7] = "106";
			PATH[8] = "110";
			PATH[9] = "111";
			PATH[10] = "112";
			PATH[11] = "113";
			PATH[12] = "114";
			PATH[13] = "115";
			PATH[14] = "116";
			break;
		case 12:
			PATH[0] = "0      ";
			PATH[1] = "10     ";
			PATH[2] = "110    ";
			PATH[3] = "1110   ";
			PATH[4] = "11110  ";
			PATH[5] = "111110 ";
			PATH[6] = "1111110";
			PATH[7] = "1111111";
			break;
		case 13:
			PATH[0] = "0000";
			PATH[1] = "0001";
			PATH[2] = "0002";
			PATH[3] = "0003";
			PATH[4] = "0004";
			PATH[5] = "0005";
			PATH[6] = "0006";
			PATH[7] = "0007";
			
			PATH[8] = "0010";
			PATH[9] = "0011";
			PATH[10] = "0012";
			PATH[11] = "0013";
			PATH[12] = "0014";
			PATH[13] = "0015";
			PATH[14] = "0016";
			PATH[15] = "0017";
			
			PATH[16] = "01 0";
			PATH[17] = "01 1";
			PATH[18] = "01 2";
			PATH[19] = "01 3";
			PATH[20] = "01 4";
			PATH[21] = "01 5";
			PATH[22] = "01 6";
			PATH[23] = "01 7";
			
			PATH[24] = "1000";
			PATH[25] = "1001";
			PATH[26] = "1002";
			PATH[27] = "1003";
			PATH[28] = "1004";
			PATH[29] = "1005";
			PATH[30] = "1006";
			PATH[31] = "1007";
			
			PATH[32] = "1010";
			PATH[33] = "1011";
			PATH[34] = "1012";
			PATH[35] = "1013";
			PATH[36] = "1014";
			PATH[37] = "1015";
			PATH[38] = "1016";
			PATH[39] = "1017";
			
			PATH[40] = "11 0";
			PATH[41] = "11 1";
			PATH[42] = "11 2";
			PATH[43] = "11 3";
			PATH[44] = "11 4";
			PATH[45] = "11 5";
			PATH[46] = "11 6";
			PATH[47] = "11 7";
			break;
		case 14:
			PATH[0] = "00";
			PATH[1] = "01";
			PATH[2] = "1 ";
			break;
		case 15:
			PATH[0] = "000    ";
			PATH[1] = "001    ";
			PATH[2] = "01     ";
			PATH[3] = "1  0   ";
			PATH[4] = "1  1000";
			PATH[5] = "1  1001";
			PATH[6] = "1  1010";
			PATH[7] = "1  1011";
			PATH[8] = "1  1100";
			PATH[9] = "1  1101";
			PATH[10] = "1  1110";
			PATH[11] = "1  1111";
			break;
		case 16:
			PATH[0] = "00";
			PATH[1] = "01";
			PATH[2] = "1 ";
			break;
		case 17:
			PATH[0] = "00";
			PATH[1] = "01";
			PATH[2] = "10";
			PATH[3] = "11";
			break;
		case 18:
			PATH[0] = "0   0";
			PATH[1] = "0   1";
			PATH[2] = "10  0";
			PATH[3] = "10  1";
			PATH[4] = "110 0";
			PATH[5] = "110 1";
			PATH[6] = "11100";
			PATH[7] = "11101";
			PATH[8] = "11110";
			PATH[9] = "11111";
			break;
		case 19:
			PATH[0] = "0011";
			PATH[1] = "0101";
			PATH[2] = "0110";
			PATH[3] = "1   ";
			break;
		case 20:
			PATH[0] = "0";
			PATH[1] = "1";
			break;
		case 21:
			PATH[0] = "0";
			PATH[1] = "1";
			break;
		case 22:
			PATH[0] = "00";
			PATH[1] = "01";
			PATH[2] = "10";
			PATH[3] = "11";
			break;
		case 23:
			PATH[0] = "0  ";
			PATH[1] = "10 ";
			PATH[2] = "110";
			PATH[3] = "111";
			break;
		case 24:
			PATH[0] = "00   ";
			PATH[1] = "0100 ";
			PATH[2] = "0101 ";
			PATH[3] = "011 0";
			PATH[4] = "011 1";
			PATH[5] = "10   ";
			PATH[6] = "1100 ";
			PATH[7] = "1101 ";
			PATH[8] = "111 0";
			PATH[9] = "111 1";
			break;
		case 25:
			PATH[0] = "0";
			PATH[1] = "1";
			PATH[2] = "2";
			PATH[3] = "3";
			PATH[4] = "4";
			PATH[5] = "5";
			PATH[6] = "6";
			PATH[7] = "7";
			PATH[8] = "8";
			PATH[9] = "9";
			PATH[10] = String.valueOf((char)58);
			PATH[11] = String.valueOf((char)59);
			PATH[12] = String.valueOf((char)60);
			PATH[13] = String.valueOf((char)61);
			PATH[14] = String.valueOf((char)62);
			break;
		case 26:
			PATH[0] = "0 ";
			PATH[1] = "10";
			PATH[2] = "11";
			break;
		case 27:
			PATH[0] = "00        ";
			PATH[1] = "01        ";
			PATH[2] = "1 00      ";
			PATH[3] = "1 01      ";
			PATH[4] = "1 1 00    ";
			PATH[5] = "1 1 01    ";
			PATH[6] = "1 1 1 00  ";
			PATH[7] = "1 1 1 01  ";
			PATH[8] = "1 1 1 1 00";
			PATH[9] = "1 1 1 1 01";
			PATH[10] = "1 1 1 1 1 ";
			break;
		case 28:
			PATH[0] = "0   ";
			PATH[1] = "10  ";
			PATH[2] = "11  ";
			PATH[3] = "12  ";
			PATH[4] = "2 00";
			PATH[5] = "2 01";
			PATH[6] = "2 02";
			PATH[7] = "2 1 ";
			PATH[8] = "2 2 ";
			PATH[9] = "2 3 ";
			PATH[10] = "3   ";
			break;
		}
	}
	
	public static void nodeiscoverage(int[] x , int func_num)  //visit[][2]:To sign whether the 'Yes' branch or the 'No' branch of each node has been covered with test case (that we have obtained).   
	{
		if(func_num == 1)
		{
			char tuple[]= new char[R];
			for(int i=0;i<R;i++)
				tuple[i] = (char) x[i];
			
			if(tuple[0]=='O'&&tuple[1]=='l'&&tuple[2]=='d')
				visit[0][0] = true;
			else
				visit[0][1] = true;
		}
		if(func_num == 2)
		{
			int entityId = (int)x[0];
			double delay = x[1];
			
			if(entityId<0)
			{
				visit[0][0] = true;
				visit[3][0] = true;
			}
			else
			{
				visit[0][1] = true;
				visit[3][1] = true;
			}
			if(delay<0)
				visit[1][0] = true;
			else
				visit[1][1] = true;
			if(entityId>=999999)
				visit[2][0] = true;
			else
				visit[2][1] = true;
			if(entityId!=1)
				visit[4][0] = true;
			else
				visit[4][1] = true;
		}
		if(func_num == 3)
		{
			int eventTime = (int)x[0];
			int type = (int)x[1];
			int dest = (int)x[2];
			int state = (int)x[3];
			int p = (int)x[4];
			int tag = (int)x[5];
			int src = (int)x[6];
			
			if(eventTime<50) visit[0][0] = true;
			else visit[0][1] = true;
			
			if(type==0) visit[1][0] = true;
			else if(type==1) visit[1][1] = true;
			else if(type==2) visit[1][2] = true;
			else visit[1][3] = true;
			
			if(dest<0) visit[2][0] = true;
			else visit[2][1] = true;
			
			if(src<0) visit[3][0] = true;
			else visit[3][1] = true;
			
			if(state==1) visit[4][0] = true;
			else visit[4][1] = true;
			
			if(p==0||tag==9999) visit[5][0] = true;
			else visit[5][1] = true;
		}
		if(func_num == 4)
		{
			int Direction = (int)x[0];
			String map1 = String.valueOf((char)x[1])+String.valueOf((char)x[2])+String.valueOf((char)x[3]);
			String map2 = String.valueOf((char)x[4])+String.valueOf((char)x[5])+String.valueOf((char)x[6]);
			
			if(Direction==1)
				visit[0][0] = true;
			else 
				visit[0][1] = true;
			
			if(map1.equals("001")||map1.equals("002")||map1.equals("003"))
				visit[1][0] = true;
			else
				visit[1][1] = true;
			
			if(map2.equals("004")||map2.equals("005")||map2.equals("006"))
				visit[2][0] = true;
			else
				visit[2][1] = true;
		}
		if(func_num == 5)
		{
			boolean isFinished;
			boolean cloudletCompleted;
			if((int)x[0]==0) isFinished=true;
			else isFinished=false;
			if((int)x[R-1]==0) cloudletCompleted=true;
			else cloudletCompleted=false;
			String cl = String.valueOf((char)x[1])+String.valueOf((char)x[2])+String.valueOf((char)x[3]);
		
			if(isFinished)
				visit[0][0] = true;
			else visit[0][1] = true;
			
			if(cl.equals("001")||cl.equals("002")||cl.equals("003"))
				visit[1][0] = true;
			else visit[1][1] = true;
			
			if(cloudletCompleted)
				visit[2][0] = true;
			else visit[2][1] = true;
		}
		if(func_num == 6)
		{
			String edge = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2]); 
			String pair = String.valueOf((char)x[3])+String.valueOf((char)x[4])+String.valueOf((char)x[5]); 
			boolean canSelect;
			if((int)x[6]==0) canSelect = true;
			else canSelect = false;
			
			if(edge.equals("mod"))
				visit[0][0] = true;
			else visit[0][1] = true;
			
			if(pair.equals("001")||pair.equals("002")||pair.equals("003"))
				visit[1][0] = true;
			else visit[1][1] = true;
			
			if(canSelect) visit[2][0] = true;
			else visit[2][1] = true;
			
			if((int)x[7]==2) visit[3][0] = true;
			else visit[3][1] = true;
		}
		if(func_num == 7){
			double r,g,b,max,min,h;
			r = x[0]; g = x[1]; b = x[2];
			max = Math.max(Math.max(r, g ), b);
			min = Math.min(Math.min(r, g ), b);
			h = max-min>0.0001&&r==max?(g-b)/(max-min):0;
			
			if(max-min>0.0001)
				visit[0][0] = true;
			else visit[0][1] = true;
			
			if(r==max) visit[1][0] = true;
			else visit[1][1] = true;
			
			if(h<0.0) visit[2][0] = true;
			else visit[2][1] = true;
			
			if(g==max) visit[3][0] = true;
			else visit[3][1] = true;
			
			if(b==max) visit[4][0] = true;
			else visit[4][1] = true;
		}
		if(func_num==8){
			double h,s,v,hue;
			int i;
			h=x[0];s=x[1];v=x[2];
			hue = s==0?0:h;
			
			if((int)s==0) visit[0][0] = true;
			else visit[0][1] = true;
			
			if((int)hue==1) visit[1][0] = true;
			else visit[1][1] = true;
			
			i = (int) ( 6*((int)hue==1?0:hue));
			
			if(i==0) visit[2][0] = true;
			else if(i==1) visit[2][1] = true;
			else if(i==2) visit[2][2] = true;
			else if(i==3) visit[2][3] = true;
			else if(i==4) visit[2][4] = true;
			else if(i==5) visit[2][5] = true;
			else visit[2][6] = true;
		}
		if(func_num==9){
			String ISBN = null; int k= (int) x[11]; 
			if(k>=ub[R-1]) k = k-1;
			for(int i=0;i<11;i++)
				ISBN += (char)x[i];
			
			if(ISBN.charAt(k)==' '||ISBN.charAt(k)=='-')
				visit[0][0] = true;
			else if((ISBN.charAt(k)-'0'>=0&&ISBN.charAt(k)-'0'<=9)||ISBN.charAt(k)=='X'||ISBN.charAt(k)=='x')
				visit[0][1] = true;
			else visit[0][2] = true;
			
			if(k<10) visit[1][0] = true;
			else visit[1][1] = true;
			
			if(k==10) visit[2][0] = true;
			else visit[2][1] = true;
			
			int temp = (ISBN.charAt(10)-'X'==0||ISBN.charAt(10)-'x'==0)?10:0;
			if(temp!=10)
				temp = (ISBN.charAt(10)-'0'>=0&&ISBN.charAt(10)-'0'<=9)?ISBN.charAt(10)-'0':0;
			
			if(checksum(ISBN)%11!=(temp)) visit[3][0] = true;
			else visit[3][1] = true;
			
			if(k>0) visit[4][0] = true;
			else visit[4][1] = true;
		}
		if(func_num==10){
			String ISSN = null; int k= (int) x[9];
			if(k>=ub[R-1]) k = k-1;
			for(int i=0;i<9;i++)
				ISSN += (char)x[i];
			
			if(ISSN.charAt(k)==' '||ISSN.charAt(k)=='-')
				visit[0][0] = true;
			else if((ISSN.charAt(k)-'0'>=0&&ISSN.charAt(k)-'0'<=9)||ISSN.charAt(k)=='X'||ISSN.charAt(k)=='x')
				visit[0][1] = true;
			else visit[0][2] = true;
			
			if(k<8) visit[1][0] = true;
			else visit[1][1] = true;
			
			if(k==8) visit[2][0] = true;
			else visit[2][1] = true;
			
			int temp = (ISSN.charAt(8)-'X'==0||ISSN.charAt(8)-'x'==0)?10:0;
			if(temp!=10)
				temp = (ISSN.charAt(8)-'0'>=0&&ISSN.charAt(8)-'0'<=9)?ISSN.charAt(10)-'0':0;
			
			if(checksum(ISSN)%11!=(temp)) visit[3][0] = true;
			else visit[3][1] = true;
			
			if(k>0) visit[4][0] = true;
			else visit[4][1] = true;
		}
		if(func_num==11){
			double h =  x[0];double w = x[1];double b = x[2];
			int i;
			h = 6.0* h/360.0;
			double temp = Math.floor(h);
			i = (int) temp;
			
			if(temp==-1.0) visit[0][0] = true;
			else visit[0][1] = true;
			
			if(temp==1.0) visit[1][0] = true;
			else visit[1][1] = true;

			if(i==0) visit[2][0] = true;
			else if(i==1) visit[2][1] = true;
			else if(i==2) visit[2][2] = true;
			else if(i==3) visit[2][3] = true;
			else if(i==4) visit[2][4] = true;
			else if(i==5) visit[2][5] = true;
			else visit[2][6] = true;
		}
		if(func_num==12){
			int nSample = (int)x[0], nPixel = (int) x[1];
			
			if(nSample==1&&nPixel==1) visit[0][0] = true;
			else visit[0][1] = true;
			
			if(nSample==1&&nPixel==2) visit[1][0] = true;
			else visit[1][1] = true;
			
			if(nSample==1&&nPixel==4) visit[2][0] = true;
			else visit[2][1] = true;
			
			if(nSample==2&&nPixel==2) visit[3][0] = true;
			else visit[3][1] = true;
			
			if(nSample==2&&nPixel==32) visit[4][0] = true;
			else visit[4][1] = true;
			
			if(nSample==3&&nPixel==4) visit[5][0] = true;
			else visit[5][1] = true;
			
			if(nSample==3&&nPixel==8) visit[6][0] = true;
			else visit[6][1] = true;
		}
		if(func_num == 13)
		{
			String option = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2]); 
			String extraOption = String.valueOf((char)x[3])+String.valueOf((char)x[4])+String.valueOf((char)x[5]); 
			int type = x[6];
			
			if(option.equals("   "))
				visit[0][0] = true;
			else visit[0][1] = true;
			
			if(!extraOption.equals("   "))
				visit[1][0] = true;
			else visit[1][1] = true;
			
			if(extraOption.endsWith(","))
				visit[2][0] = true;
			else visit[2][1] = true;
			
			if(type == 0) visit[3][0] = true;
			else if(type ==1) visit[3][1] = true;
			else if(type ==2) visit[3][2] = true;
			else if(type ==3) visit[3][3] = true;
			else if(type ==4) visit[3][4] = true;
			else if(type ==5) visit[3][5] = true;
			else if(type ==6) visit[3][6] = true;
			else if(type >6) visit[3][7] = true;
		}
		if(func_num == 14){
			String xmlTags = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2]); 
			String sentence = String.valueOf((char)x[3])+String.valueOf((char)x[4])+String.valueOf((char)x[5]); 
			
			if(xmlTags.equals("   "))
				visit[0][0] = true;
			else visit[0][1] = true;
			
			if(sentence.equals("   "))
				visit[1][0] = true;
			else visit[1][1] = true;
		}
		if(func_num == 15){
			Boolean nlSplitting,whitespaceTokenization;
			if(x[0]==0) nlSplitting = true;
			else nlSplitting = false;
			if(x[1]==0) whitespaceTokenization = true;
			else whitespaceTokenization = false;
			String line = String.valueOf((char)x[2])+String.valueOf((char)x[3]); 
			String isOneSentence = String.valueOf((char)x[4])+String.valueOf((char)x[5])+String.valueOf((char)x[6])+String.valueOf((char)x[7]); 
			char token,bound1,bound2;
			token = (char)x[8]; bound1 = (char)x[9]; bound2 = (char)x[10];
			
			if(nlSplitting) visit[0][0] = true;
			else visit[0][1] = true;
			
			if(whitespaceTokenization) visit[1][0] = true;
			else visit[1][1] = true;
			
			if(line.equals("/n")) visit[2][0] = true;
			else visit[2][1] = true;
			
			if(isOneSentence.equals("true")) visit[3][0] = true;
			else visit[3][1] = true;
			
			if(token==' ') visit[4][0] = true;
			else visit[4][1] = true;
			
			if(bound1==' ') visit[5][0] = true;
			else visit[5][1] = true;
			
			if(bound2==' ') visit[6][0] = true;
			else visit[6][1] = true;
		}
		if(func_num == 16){
			String annotation = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2]); 
			int nThreads = x[3];
			
			if(annotation.equals("001")||annotation.equals("002")||annotation.equals("003"))
				visit[0][0] = true;
			else
				visit[0][1] = true;
			
			if(nThreads == 1)
				visit[1][0] = true;
			else 
				visit[1][1] = true;
		}
		if(func_num == 17){
			String nerLanguage = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2])+
					String.valueOf((char)x[3])+String.valueOf((char)x[4])+String.valueOf((char)x[5])+String.valueOf((char)x[6]); 
			Boolean augment;
			if(x[7]==0) augment = true;
			else augment = false;
			
			if(nerLanguage.equals("CHINESE"))
				visit[0][0] = true;
			else visit[0][1] = true;
			
			if(augment) visit[1][0] = true;
			else visit[1][1] = true;
		}
		if(func_num==18){
			String trueCase = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2])+String.valueOf((char)x[3])+String.valueOf((char)x[4]);
			boolean overwriteText;
			if(x[5]==0) overwriteText = true;
			else overwriteText = false;
			
			if(trueCase.equals("UPPER")) visit[0][0] = true;
			else visit[0][1] = true;
			
			if(trueCase.equals("LOWER")) visit[1][0] = true;
			else visit[1][1] = true;
			
			if(trueCase.equals("INIT_")) visit[2][0] = true;
			else visit[2][1] = true;
			
			if(trueCase.equals("O    ")) visit[3][0] = true;
			else visit[3][1] = true;
			
			if(overwriteText) visit[4][0] = true;
			else visit[4][1] = true;
		}
		if(func_num == 19)   //Triangle
		{
			int a = x[0] ;
			int b = x[1] ;
			int c = x[2] ;		
					
			if((a<(b+c)) && (b<(a+c)) && (c<(a+b)))
			{	
				visit[0][0] = true ;
    			if ( ((a==b) && (a!=c))	 || ((a==c)&&(a!=b)) || ((b==c)&&(b!=a)) )	  
    				visit[1][0] = true ;
    			else
    				visit[1][1] = true ;
    			if ( (a==b) && (a==c) )
    				visit[2][0] = true  ;
    			else
    				visit[2][1] = true ;
    			if ( (a!=b) && (a!=c) && (b!=c) )
    				visit[3][0] = true ;
    			else 
    				visit[3][1] = true ;
    		}
			else
				visit[0][1] = true ;   		
		}
		
		if(func_num == 20) //Factorial
		{
			int a = x[0] ;
			
			if(a==1)
				visit[0][0] = true ;
			else
				visit[0][1] = true ;
		}
		
		if(func_num == 21)  //bubble sorting
		{
			int i1,j1;
			int[] a = new int[R];
		   	  for(i1=0;i1<R;i1++)
		   	       a[i1] = x[i1];
		   	for(j1=0;j1<=R-1;j1++) 
		   	{
		   		 for (i1=0;i1<R-1-j1;i1++)
		   		 {
		   			 if(a[i1]>a[i1+1])
		   			   { visit[0][0] = true ; break ;}
		   		 }
		   		 if(visit[0][0]) break;	   		  
		   	}
		   	if(!visit[0][0])
		   		visit[0][1] = true ;
		}
		
		if(func_num == 22)   //GCD
		{
			int m = x[0] ;
			int n = x[1] ;
			
			if (m<n)
	   		{
				visit[0][0] = true ;
	   			int t = n;
	   			n = m;
	   			m = t;
	   		}
			else
				visit[0][1] = true ;
			
			int r;
		   	r = m % n;
		    m = n;
		   	n = r;
		   	
		   	while (r != 0)
	   		{
		   		visit[1][0] = true ;
	   		    r = m % n;
	   		    m = n;
	   		    n = r;
	   		    break ;
	   	    }
		   	
		   	if(!visit[1][0])
		   		visit[1][1] = true ;
		}
		
		if(func_num == 23)  //Middle
		{
			int a = x[0] ;
			int b = x[1] ;
			int c = x[1] ;
			
			if( ( (a < b) && (b < c) ) || ((c<b) && (b<a)) )
	   			visit[0][0] = true ;
	   		else if ( ( (a < c) && (c < b) ) || ((b<c) && (c<a)) )
	   			{visit[1][0] = true ; visit[0][1] = true ;}
	   		else if ( ( (b < a) && (a < c) ) || ((c<a) && (a<b)) )
	   			{visit[2][0] = true ; visit[0][1] = true ; visit[1][1] = true ;}
	   		else
	   			{visit[0][1] = true ; visit[1][1] = true ; visit[2][1] = true ;}
		}
		
		if(func_num == 24)  //Tomorrow
		{
			int Day = x[0] ;
      	    int Year = x[1] ;
      	    int Month = x[2] ;
      	    int Date = x[3] ;
      	    
      	    if (Day == 7)
	   	    	visit[0][0] = true ;
	   	    else
	   	    	visit[0][1] = true ;
      	  
	      	if (Month == 12 && Date == 31)
	   	    {
	   	    	visit[1][0] = true ;
	   	    }		   	    		   	    
	   	    else if(Month == 2 && Date == 28)
	   	    {	
	   	    	visit[1][1] = true ;
	   	    	visit[2][0] = true ;
	   	    	if(isRun(Year))
	   	    		visit[3][0] = true ;
	   	    	else
	   	    		visit[3][1] = true ;
	   	    }		   	    
	   	    else if((Month != 12  && Date == 31) || (Month == 2 && Date == 29) 
	   	    		|| ((Month == 4 || Month == 6 || Month == 9 || Month == 11) && Date == 30))
	   	    {
	   	    	visit[1][1] = true ;
	   	    	visit[2][1] = true ;
	   	    	visit[4][0] = true ;
	   	    }
	   	    else 
	   	    {
	   	    	visit[1][1] = true ;
	   	    	visit[2][1] = true ;
	   	    	visit[4][1] = true ;
	   	    }
		}
		if(func_num == 25) //calculator
	    {
	    	 int ch2 = x[1] ;
	    	 char cmd = (char) (ch2) ;
			
			 switch(cmd)
	    	 {
	    	 case '+' : visit[0][0] = true ; break ;
	    	 case '-' : visit[0][1] = true ; break ;
	    	 case '.' : visit[0][2] = true ; break ;
	    	 case '/' : visit[0][3] = true ; break ;
	    	 case 'p' : visit[0][4] = true ; break ;
	    	 case 'a' : visit[0][5] = true ; break ;
	    	 case 'b' : visit[0][6] = true ; break ;
	    	 case 'c' : visit[0][7] = true ; break ;
	    	 case 'd' : visit[0][8] = true ; break ;
	    	 case 'e' : visit[0][9] = true ; break ;
	    	 case 'f' : visit[0][10] = true ; break ;
	    	 case 'S' : visit[0][11] = true ; break ;
	    	 case 'C' : visit[0][12] = true ; break ;
	    	 case 'T' : visit[0][13] = true ; break ;
	    	 default  : visit[0][14] = true ; break ;
	    	 }
		}
		if(func_num == 26) //commission
		{
			int totallocks = x[0] ;
			int totalstocks = x[1] ;
			int totalbarrels = x[2] ;
			
			double  lockprice = 45.0 ;
			double  stockprice = 30.0 ;
			double  barrelprice = 25.0 ;
			
			double  locksales = lockprice * totallocks ;
			double  stocksales = stockprice * totalstocks ;
			double  barrelsales = barrelprice * totalbarrels ;
			double  sales = locksales + stocksales + barrelsales ;

			if(sales > 1800.0)
			  visit[0][0] = true;
			else if(sales > 500.0)
			{
				visit[0][1] = true ;
				visit[1][0] = true ;
			}
			else
			{
			  visit[0][1] = true ;
			  visit[1][1] = true ;
			}
		}
		if(func_num == 27) //premium
		{
			int  driverage = x[0] ;
			int  points = x[1] ;
			
			if(driverage >=16 && driverage < 20)
			{
			    visit[0][0] = true ;
			    if(points <= 1)
			      visit[1][0] = true ;
			    else
			    	visit[1][1] = true ;
			}
			else if(driverage >= 20 && driverage < 25)
			{
				visit[0][1] = true ;
				visit[2][0] = true ;
			    if(points < 3)
			       visit[3][0] = true ; 
			    else
			       visit[3][1] = true ; 
			}
			else if(driverage >= 25 && driverage < 45)
			{
				visit[0][1] = true ;
				visit[2][1] = true ;
				visit[4][0] = true ;
			    if(points < 5)
			       visit[5][0] = true ; 
			    else
			       visit[5][1] = true ; 
			}
			else if(driverage >= 45 && driverage < 60)
			{
				visit[0][1] = true ;
				visit[2][1] = true ;
				visit[4][1] = true ;
				visit[6][0] = true ;
			    if(points < 7)
			    	visit[7][0] = true  ;
			    else
			    	visit[7][1] = true ;
			}
			else if(driverage >= 60 && driverage < 100)
			{
				visit[0][1] = true ;
				visit[2][1] = true ;
				visit[4][1] = true ;
				visit[6][1] = true ;
				visit[8][0] = true ;
			    if(points < 5)
			    	visit[9][0] = true ;
			    else
			    	visit[9][1] = true ;
			}
			else
			{
				visit[0][1] = true ;
				visit[2][1] = true ;
				visit[4][1] = true ;
				visit[6][1] = true ;
				visit[8][1] = true ;
			}
		}
		
		if(func_num == 28) //decision
		{
			int a = x[0] ;
			int b = x[1] ;
			int c = x[2] ;	
			int d = x[3] ;
					
			if(a==3971)
				visit[0][0] = true ;
			else if(a==5085)
			{
				visit[0][1] = true ;
				if(c==5448)
					visit[1][0] = true ;
				else if(c==2463)
					visit[1][1] = true ;
				else
					visit[1][2] = true ;
			}
			else if(a==5174)
			{
				visit[0][2] = true ;
				if(b==4040)
				{
					visit[2][0] = true ;
					if(d==5148)
						visit[3][0] = true ;
					else if(d==4662)
						visit[3][1] = true ;
					else
						visit[3][2] = true ;
				}
				else if(b==5448)
					visit[2][1] = true ;
				else if(b==3268)
					visit[2][2] = true ;
				else
					visit[2][3] = true ;
			}
			else 
				visit[0][3] = true ;  
		}
	}
	
	public static int pathnum(int[] x , int func_num)
	{
		int path = -1;
		
		if(func_num == 1)
		{
			char tuple[]= new char[R];
			for(int i=0;i<R;i++)
				tuple[i] = (char) x[i];
			
			if(tuple[0]=='O'&&tuple[1]=='l'&&tuple[2]=='d')
				path = 0;
			else
				path = 1;
		}
		if(func_num == 2)
		{
			int entityId = (int)x[0];
			double delay = x[1];
			
			if(entityId<0)
			{
				path = 0;
			}else{
				if(delay<0)
				{
					delay = 0;
					if(delay>=999999)
						path = 1;
					else if(entityId<0)
						path = 2;
					else if(entityId!=1)
						path = 3;
					else 
						path = 4;
						
				}
				else{
					if(delay>=999999)
						path = 5;
					else if(entityId<0)
						path = 6;
					else if(entityId!=1)
						path = 7;
					else 
						path = 8;
				}
			}
		}
		if(func_num == 3)
		{
			int eventTime = (int) x[0];
			int type = (int)x[1];
			int dest = (int)x[2];
			int state = (int)x[3];
			int p = (int)x[4];
			int tag = (int)x[5];
			int src = (int)x[6];
			
			if(eventTime<50)
				path = 0;
			else{
				if(type==0)
					path = 1;
				else if(type==3)
					path = 8;
				else if(type==2){
					if(src<0)
						path = 6;
					else
						path = 7;
				}else{
					if(dest<0)
						path = 2;
					else if(state==1){
						path = 5;
					}else{
						if(p==0||tag==9999)
							path = 3;
						else
							path = 4;
					}
				}
			}
		}
		if(func_num == 4)
		{
			int Direction = (int)x[0];
			String map1 = String.valueOf((char)x[1])+String.valueOf((char)x[2])+String.valueOf((char)x[3]);
			String map2 = String.valueOf((char)x[4])+String.valueOf((char)x[5])+String.valueOf((char)x[6]);
			
			if(Direction==1){
				if(map1.equals("001")||map1.equals("002")||map1.equals("003")){
					if(map2.equals("004")||map2.equals("005")||map2.equals("006"))
						path = 0;
					else
						path = 1;
				}else{
					if(map2.equals("004")||map2.equals("005")||map2.equals("006"))
						path = 2;
					else
						path = 3;
				}
			}
			else path = 4;
		}
		if(func_num == 5)
		{
			boolean isFinished;
			boolean cloudletCompleted;
			if((int)x[0]==0) isFinished=true;
			else isFinished=false;
			if((int)x[R-1]==0) cloudletCompleted=true;
			else cloudletCompleted=false;
			String cl = String.valueOf((char)x[1])+String.valueOf((char)x[2])+String.valueOf((char)x[3]);
			
			if(isFinished)
			{
				if(cl.equals("001")||cl.equals("002")||cl.equals("003"))
					if(cloudletCompleted)
						path = 0;
					else
						path = 1;
				else
					if(cloudletCompleted)
						path = 2;
					else
						path = 3;
			}else
				if(cloudletCompleted)
					path = 4;
				else
					path = 5;
		}
		if(func_num == 6)
		{
			String edge = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2]); 
			String pair = String.valueOf((char)x[3])+String.valueOf((char)x[4])+String.valueOf((char)x[5]); 
			boolean canSelect;
			if((int)x[6]==0) canSelect = true;
			else canSelect = false;
			
			if(edge.equals("mod"))
			{
				if(pair.equals("001")||pair.equals("002")||pair.equals("003")){
					if(canSelect){
						if((int)x[7]==2)
							path = 0;
						else
							path = 1;
					}else
						path = 2;
				}else{
					if(canSelect){
						if((int)x[7]==2)
							path = 3;
						else
							path = 4;
					}else
						path = 5;
				}			
			}else
				path = 6;
		}

		if(func_num == 7){
			double r,g,b,max,min; double h;
			r = x[0]; g = x[1]; b = x[2];
			max = Math.max(Math.max(r, g ), b);
			min = Math.min(Math.min(r, g ), b);
			h = max-min>0.0001&&r==max?(g-b)/(max-min):0;
			
			if(max-min>0.0001){
				if(r==max){
					if(h<0.0)
						path = 0;
					else
						path = 1;
				}else if(g==max)
					path = 2;
				else if(b==max)
					path = 3;
				else 
					path = 4;
			}else
				path = 5;
		}
		if(func_num==8){
			double h,s,v,hue;
			int i;
			h=x[0];s=x[1];v=x[2];
			hue = s==0?0:h;
			
			if((int)s==0) path = 0;
			else{
				hue = h;
				if((int)hue==1){
					hue =0.0;
					i = (int) hue;
					if(i==0)
						path = 1;
					else if(i==1)
						path = 2;
					else if(i==2)
						path = 3;
					else if(i==3)
						path = 4;
					else if(i==4)
						path = 5;
					else if(i==5)
						path = 6;
					else 
						path = 7;
				}else{
					hue *=6.0;
					i = (int) ( 6*((int)hue==1?0:hue));
					if(i==0)
						path = 8;
					else if(i==1)
						path = 9;
					else if(i==2)
						path = 10;
					else if(i==3)
						path = 11;
					else if(i==4)
						path = 12;
					else if(i==5)
						path = 13;
					else 
						path = 14;
				}
			}
		}
		if(func_num==9){
			String ISBN = null; int k= (int) x[R-1];
			
			for(int i=0;i<R-1;i++)
				ISBN += (char)x[i];
			
//			System.out.println(k);
			
			if(ISBN.charAt(k)==' '||ISBN.charAt(k)=='-')
				path = 0;
			else if((ISBN.charAt(k)-'0'>=0&&ISBN.charAt(k)-'0'<=9)||ISBN.charAt(k)=='X'||ISBN.charAt(k)=='x'){
				if(k<10)
					path = 1;
				else
					if(k==10){
						int temp = (ISBN.charAt(10)-'X'==0||ISBN.charAt(10)-'x'==0)?10:0;
						if(temp!=10)
							temp = (ISBN.charAt(10)-'0'>=0&&ISBN.charAt(10)-'0'<=9)?ISBN.charAt(10)-'0':0;
						if(checksum(ISBN)%11!=temp)
							path = 2;
						else
							path = 3;
					}else
						path=4;
			}
			else{
				if(k>0)
					path = 5;
				else 
					path = 6;
			}
		}
		if(func_num==10){
			String ISSN = null; int k= (int) x[9];
			if(k>=ub[R-1]) k = k-1;
			for(int i=0;i<9;i++)
				ISSN += (char)x[i];
			
			if(ISSN.charAt(k)==' '||ISSN.charAt(k)=='-')
				path = 0;
			else if((ISSN.charAt(k)-'0'>=0&&ISSN.charAt(k)-'0'<=9)||ISSN.charAt(k)=='X'||ISSN.charAt(k)=='x'){
				if(k<8)
					path = 1;
				else
					if(k==8){
						int temp = (ISSN.charAt(8)-'X'==0||ISSN.charAt(8)-'x'==0)?10:0;
						if(temp!=10)
							temp = (ISSN.charAt(8)-'0'>=0&&ISSN.charAt(8)-'0'<=9)?ISSN.charAt(8)-'0':0;
						if(checksum(ISSN)%11!=temp)
							path = 2;
						else
							path = 3;
					}else
						path=4;
			}
			else{
				if(k>0)
					path = 5;
				else 
					path = 6;
			}
		}
		if(func_num==11){
			double h =  x[0];double w = x[1];double b = x[2];
			int i;
			h = 6.0* h/360.0;
			double temp = Math.floor(h);
			i = (int) temp;
			
			if(temp==-1.0)
				path = 0;
			else{
				if(temp==1.0){
					if(i==0)
						path = 1;
					else if(i==1)
						path = 2;
					else if(i==2)
						path = 3;
					else if(i==3)
						path = 4;
					else if(i==4)
						path = 5;
					else if(i==5)
						path = 6;
					else
						path = 7;
				}else{
					if(i==0)
						path = 8;
					else if(i==1)
						path = 9;
					else if(i==2)
						path = 10;
					else if(i==3)
						path = 11;
					else if(i==4)
						path = 12;
					else if(i==5)
						path = 13;
					else
						path = 14;
				}
			}
		}
		if(func_num==12){
			int nSample = (int)x[0], nPixel = (int) x[1];
			
			if(nSample==1&&nPixel==1)
				path = 0;
			else if(nSample==1&&nPixel==2)
				path = 1;
			else if(nSample==1&&nPixel==4)
				path = 2;
			else if(nSample==2&&nPixel==2)
				path = 3;
			else if(nSample==2&&nPixel==32)
				path = 4;
			else if(nSample==3&&nPixel==4)
				path = 5;
			else if(nSample==3&&nPixel==8)
				path = 6;
			else
				path = 7;
		}
		if(func_num == 13)
		{
			String option = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2]); 
			String extraOption = String.valueOf((char)x[3])+String.valueOf((char)x[4])+String.valueOf((char)x[5]); 
			int type = x[6];
			
			if(option.equals("   "))
			{
				if(!extraOption.equals("   ")){
					if(extraOption.endsWith(","))
					{//"0000"-"0007"
						if(type == 0) path = 0;
						else if(type == 1) path = 1;
						else if(type == 2) path = 2;
						else if(type == 3) path = 3;
						else if(type == 4) path = 4;
						else if(type == 5) path = 5;
						else if(type == 6) path = 6;
						else if(type >6) path = 7;
					}
					else{//"0010"-"0017"
						if(type == 0) path = 8;
						else if(type == 1) path = 9;
						else if(type == 2) path = 10;
						else if(type == 3) path = 11;
						else if(type == 4) path = 12;
						else if(type == 5) path = 13;
						else if(type == 6) path = 14;
						else if(type >6) path = 15;
					}
				}else{//"01 0"-"01 7"
					if(type == 0) path = 16;
					else if(type == 1) path = 17;
					else if(type == 2) path = 18;
					else if(type == 3) path = 19;
					else if(type == 4) path = 20;
					else if(type == 5) path = 21;
					else if(type == 6) path = 22;
					else if(type >6) path = 23;
				}
			}else{
				if(!extraOption.equals("   ")){
					if(extraOption.endsWith(","))
					{//"1000"-"1007"
						if(type == 0) path = 24;
						else if(type == 1) path = 25;
						else if(type == 2) path = 26;
						else if(type == 3) path = 27;
						else if(type == 4) path = 28;
						else if(type == 5) path = 29;
						else if(type == 6) path = 30;
						else if(type >6) path = 31;
					}
					else{//"1010"-"1017"
						if(type == 0) path = 32;
						else if(type == 1) path = 33;
						else if(type == 2) path = 34;
						else if(type == 3) path = 35;
						else if(type == 4) path = 36;
						else if(type == 5) path = 37;
						else if(type == 6) path = 38;
						else if(type >6) path = 39;
					}
				}else{//"11 0"-"11 7"
					if(type == 0) path = 40;
					else if(type == 1) path = 41;
					else if(type == 2) path = 42;
					else if(type == 3) path = 43;
					else if(type == 4) path = 44;
					else if(type == 5) path = 45;
					else if(type == 6) path = 46;
					else if(type >6) path = 47;
				}
			}

		}
		if(func_num == 14){
			String xmlTags = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2]); 
			String sentence = String.valueOf((char)x[3])+String.valueOf((char)x[4])+String.valueOf((char)x[5]); 
			
			if(xmlTags.equals("   "))
			{
				if(sentence.equals("   "))
					path = 0;
				else
					path = 1;
			}else
				path = 2;
		}
		if(func_num == 15){
			Boolean nlSplitting,whitespaceTokenization;
			if(x[0]==0) nlSplitting = true;
			else nlSplitting = false;
			if(x[1]==0) whitespaceTokenization = true;
			else whitespaceTokenization = false;
			String line = String.valueOf((char)x[2])+String.valueOf((char)x[3]); 
			String isOneSentence = String.valueOf((char)x[4])+String.valueOf((char)x[5])+String.valueOf((char)x[6])+String.valueOf((char)x[7]); 
			char token,bound1,bound2;
			token = (char)x[8]; bound1 = (char)x[9]; bound2 = (char)x[10];
			
			if(nlSplitting){
				if(whitespaceTokenization){
					if((line.equals("/n")))
						path = 0;
					else
						path = 1;
				}
				else
					path = 2;
			}
			else{
				if(isOneSentence.equals("true"))
					path = 3;
				else{
					if(token==' '){
						if(bound1==' '){
							if(bound2==' ')
								path = 4;
							else
								path = 5;
						}else{
							if(bound2==' ')
								path = 6;
							else
								path = 7;
						}
					}else{
						if(bound1==' '){
							if(bound2==' ')
								path = 8;
							else
								path = 9;
						}else{
							if(bound2==' ')
								path = 10;
							else
								path = 11;
						}
					}
				}
			}
		}
		if(func_num == 16)
		{
			String annotation = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2]); 
			int nThreads = x[3];
			
			if(annotation.equals("001")||annotation.equals("002")||annotation.equals("003"))
			{
				if(nThreads == 1)
					path =0;
				else 
					path = 1;
			}else
				path =2;
		}
		if(func_num == 17){
			String nerLanguage = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2])+
					String.valueOf((char)x[3])+String.valueOf((char)x[4])+String.valueOf((char)x[5])+String.valueOf((char)x[6]); 
			Boolean augment;
			if(x[7]==0) augment = true;
			else augment = false;
			
			if(nerLanguage.equals("CHINESE")){
				if(augment)
					path = 0;
				else 
					path = 1;
			}else{
				if(augment)
					path = 2;
				else 
					path = 3;
			}
		}
		if(func_num==18){
			String trueCase = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2])+String.valueOf((char)x[3])+String.valueOf((char)x[4]);
			boolean overwriteText;
			if(x[5]==0) overwriteText = true;
			else overwriteText = false;
			
			if(trueCase.equals("UPPER")){
				if(overwriteText)
					path = 0;
				else
					path = 1;
			}else if(trueCase.equals("LOWER")){
				if(overwriteText)
					path = 2;
				else
					path = 3;
			}else if(trueCase.equals("INIT_")){
				if(overwriteText)
					path = 4;
				else
					path = 5;
			}else if(trueCase.equals("O    ")){
				if(overwriteText)
					path = 6;
				else
					path = 7;
			}else{
				if(overwriteText)
					path = 8;
				else
					path = 9;
			}
		}
		if(func_num == 19) //Triangle
		{
			int a = x[0] ;
			int b = x[1] ;
			int c = x[2] ;		
					
			if((a<(b+c)) && (b<(a+c)) && (c<(a+b)))
			{	
    			if ( ((a==b) && (a!=c))	 || ((a==c)&&(a!=b)) || ((b==c)&&(b!=a)) )	  
    				path = 0 ;
    			if ( (a==b) && (a==c) )
    				path = 1  ;
    			if ( (a!=b) && (a!=c) && (b!=c) )
    				path = 2 ;
    		}
			else
				path = 3 ;   
		}
		
		if(func_num == 20) //Factorial
		{
			int a = x[0] ;
			
			if(a==1)
				path = 0 ;
			else
				path = 1 ;
		}
		
		if(func_num == 21)  //bubble sorting
		{
			int i1,j1;
			int[] a = new int[R];
		   	for(i1=0;i1<R;i1++)
		   	     a[i1] = x[i1];
		   	for(j1=0;j1<=R-1;j1++) 
		   	{
		   		 for (i1=0;i1<R-1-j1;i1++)
		   		 {
		   			 if(a[i1]>a[i1+1])
		   			   { path = 0 ; break ;}
		   		 }
		   		 if(path == 0) break;	   		  
		   	}
		   	if(path != 0)
		   		path = 1 ;
		}
		
		if(func_num == 22)   //GCD
		{
			int m = x[0] ;
			int n = x[1] ;
			boolean d[] =new boolean[2] ;
			
			if (m<n)
	   		{
				d[0] = true ;
	   			int t = n;
	   			n = m;
	   			m = t;
	   		}
			
			int r;
		   	r = m % n;
		    m = n;
		   	n = r;
		   	
		   	while (r != 0)
	   		{
		   		d[1] = true ;
	   		    r = m % n;
	   		    m = n;
	   		    n = r;
	   		    break ;
	   	    }
		   	
		   	if(d[0] && d[1])
		   		path = 0 ;
		   	else if(d[0] && (!d[1]))
		   		path = 1 ;
		   	else if((!d[0]) && d[1])
		   		path = 2 ;
		   	else
		   		path =3 ;
		}
		
		if(func_num == 23)  //Middle
		{
			  int a = x[0] ;
		      int b = x[1] ;
		   	  int c = x[2] ;
		   	  		  		
	   		  if( ( (a < b) && (b < c) ) || ((c<b) && (b<a)) )
	   			  path = 0 ;	   		    
	   		  else if ( ( (a < c) && (c < b) ) || ((b<c) && (c<a)) )
	   			  path = 1 ;
	   		  else if ( ( (b < a) && (a < c) ) || ((c<a) && (a<b)) )
	   			  path = 2 ;
	   		  else
	   			  path = 3 ;
		}
		
		if(func_num == 24)  //Tomorrow
		{
			int Day = x[0] ;
      	    int Year = x[1] ;
      	    int Month = x[2] ;
      	    int Date = x[3] ;
      	    
      	    if (Day == 7)
      	    {
      	     	if (Month == 12 && Date == 31)
      	     		path = 0 ;
    	   	    else if(Month == 2 && Date == 28)
    	   	    {
    	   	    	if(isRun(Year))
    	   	    		path = 1 ;
    	   	    	else
    	   	    		path = 2 ;
    	   	    }
    	   	    	
    	   	    else if((Month != 12  && Date == 31) || (Month == 2 && Date == 29) 
    	   	    		|| ((Month == 4 || Month == 6 || Month == 9 || Month == 11) && Date == 30))
    	   	    	path = 3 ;
    	   	    else 
    	   	    	path = 4 ;
      	    }
      	    else
      	    {
      	    	if (Month == 12 && Date == 31)
      	     		path = 5 ;
    	   	    else if(Month == 2 && Date == 28)
    	   	    {
    	   	    	if(isRun(Year))
    	   	    		path = 6 ;
    	   	    	else
    	   	    		path = 7 ;
    	   	    }
    	   	    	
    	   	    else if((Month != 12  && Date == 31) || (Month == 2 && Date == 29) 
    	   	    		|| ((Month == 4 || Month == 6 || Month == 9 || Month == 11) && Date == 30))
    	   	    	path = 8 ;
    	   	    else 
    	   	    	path = 9 ;
      	    }     	  
		}
		if(func_num == 25) //calculator
	     {
	    	 int ch2 = x[1] ;
	    	 char cmd = (char) (ch2) ;
	    	 
	    	 switch(cmd)
	    	 {
	    	 case '+' : path = 0 ; break ;
	    	 case '-' : path = 1 ; break ;
	    	 case '.' : path = 2 ; break ;
	    	 case '/' : path = 3 ; break ;
	    	 case 'p' : path = 4 ; break ;
	    	 case 'a' : path = 5 ; break ;
	    	 case 'b' : path = 6 ; break ;
	    	 case 'c' : path = 7 ; break ;
	    	 case 'd' : path = 8 ; break ;
	    	 case 'e' : path = 9 ; break ;
	    	 case 'f' : path = 10 ; break ;
	    	 case 'S' : path = 11 ; break ;
	    	 case 'C' : path = 12 ; break ;
	    	 case 'T' : path = 13 ; break ;
	    	 default : path = 14 ; break ;
	    	 }
	     }
		if(func_num == 26) //commission
		{
			int totallocks = x[0] ;
			int totalstocks = x[1] ;
			int totalbarrels = x[2] ;
			
			double  lockprice = 45.0 ;
			double  stockprice = 30.0 ;
			double  barrelprice = 25.0 ;
			
			double  locksales = lockprice * totallocks ;
			double  stocksales = stockprice * totalstocks ;
			double  barrelsales = barrelprice * totalbarrels ;
			double  sales = locksales + stocksales + barrelsales ;

			if(sales > 1800.0)
			  path = 0;
			else if(sales > 500.0)
			  path = 1 ;
			else
	          path = 2 ;
		}
		if(func_num == 27) //premium
		{
			int  driverage = x[0] ;
			int  points = x[1] ;
			
			if(driverage >=16 && driverage < 20)
			{			 
			    if(points <= 1)
			      path = 0 ;
			    else
			      path = 1 ;
			}
			else if(driverage >= 20 && driverage < 25)
			{
			    if(points < 3)
			       path = 2 ; 
			    else
			       path = 3 ; 
			}
			else if(driverage >= 25 && driverage < 45)
			{
			    if(points < 5)
			       path = 4 ; 
			    else
			       path = 5 ; 
			}
			else if(driverage >= 45 && driverage < 60)
			{
			    if(points < 7)
			    	path = 6  ;
			    else 
			    	path = 7 ;
			}
			else if(driverage >= 60 && driverage < 100)
			{
			    if(points < 5)
			    	path = 8 ;
			    else
			    	path = 9 ;
			}
			else
                path = 10 ;
		}
		
		if(func_num == 28) //decision
		{
			int a = x[0] ;
			int b = x[1] ;
			int c = x[2] ;	
			int d = x[3] ;
					
			if(a==3971)
				path=0;
			else if(a==5085)
			{
				if(c==5448)
					path=1;
				else if(c==2463)
					path=2;
				else
					path=3;
			}
			else if(a==5174)
			{
				if(b==4040)
				{
					if(d==5148)
						path=4;
					else if(d==4662)
						path=5;
					else
						path=6;
				}
				else if(b==5448)
					path=7;
				else if(b==3268)
					path=8;
				else
				    path=9;
			}
			else 
				path = 10 ;   
		}
		
		return path;
	}
	
	public static double benchmarkfunction (int[] x , int func_num, int path_num)
	{
		double[] fit = new double[NODENUM] ;  //fit[k]表示测试用例经过节点k时，在该节点的适应值
		double[][] f = new double[NODENUM][BRANCH];
		double Fitness = 0 ;    //测试用例的适应值
		
		if(func_num == 1)
		{
			char tuple[]= new char[R];
			for(int i=0;i<R;i++)
				tuple[i] = (char) x[i];
			
			double v1=0;
			
			if(tuple[0]=='O'&&tuple[1]=='l'&&tuple[2]=='d')
				v1 = 0;
			else
				v1 = Math.abs(tuple[0]-'O')+K+ Math.abs(tuple[1]-'l')+K + Math.abs(tuple[2]-'d')+K;
			f[0][0] = v1;
			
			if(tuple[0]!='O'||tuple[1]!='l'||tuple[2]!='d')
				v1 = 0;
			else{
				v1 = Math.min(Math.abs(tuple[0]-'O'), Math.abs(tuple[1]-'l'));
				v1 = Math.min(v1, Math.abs(tuple[2]-'d'));
			}
			f[0][1] = v1;
		}
		if(func_num == 2)
		{
			int entityId = (int)x[0];
			double delay = x[1];
			
			double v1=0,v2=0;
			
			//分支节点1
			if(entityId<0) v1=0;
			else v1 = entityId+K;
			f[0][0] = v1;
			
			if(entityId>=0) v2=0;
			else v2 = 0-entityId+K;
			f[0][1] = v2;
			
			//分支节点2
			if(delay<0) {v1=0;delay=0;}
			else v1 = delay+K;
			f[1][0] = v1;
			
			if(delay>=0) v2=0;
			else v2 = -delay+K;
			f[1][1] = v2;
				
			//分支节点3
			if(delay>=999999) v1=0;
			else v1 = 999999-delay+K;
			f[2][0] = v1;
			
			if(delay<999999) v2=0;
			else v2 = delay-999999+K;
			f[2][1] = v2;
			
			//分支节点4
			if(entityId<0) v1=0;
			else v1 = entityId+K;
			f[3][0] = v1;
			
			if(entityId>=0) v2 = 0;
			else v2 = 0-entityId+K;
			f[3][1] = v2;
			
			//分支节点5
			if(entityId!=1) v1=0;
			else v1 = K;
			f[4][0] = v1;
			
			if(entityId==1) v2 = 0;
			else v2 = Math.abs(entityId+1)+K;
			f[4][1] = v2;
		}
		if(func_num == 3)
		{
			int eventTime = (int)x[0];
			int type = (int)x[1];
			int dest = (int)x[2];
			int state = (int)x[3];
			int p = (int)x[4];
			int tag = (int)x[5];
			int src = (int)x[6];
			
			double v1=0;
			//分支节点0
			if(eventTime<50) v1=0;
			else v1 = eventTime-50+K;
			f[0][0] = v1;
			
			if(eventTime>=50) v1=0;
			else v1 = 50-eventTime+K;
			f[0][1] = v1;
			
			//分支节点1
			if(type==0)	f[1][0] = 0;
			else	f[1][0] = Math.abs(type-0)+K;
			if(type==1)	f[1][1] = 0;
			else	f[1][1] = Math.abs(type-1)+K;
			if(type==2)	f[1][2] = 0;
			else	f[1][2] = Math.abs(type-2)+K;
			if(type==3)	f[1][3] = 0;
			else	f[1][3] = Math.abs(type-3)+K;

			
			//分支节点2
			if(dest<0) v1=0;
			else v1 = dest+K;
			f[2][0] = v1;
			
			if(dest>=0) v1=0;
			else v1=-dest+K;
			f[2][1] = v1;
			
			//分支节点3
			if(src<0) v1=0;
			else v1=src+K;
			f[3][0] = v1;
			
			if(src>=0) v1=0;
			else v1= -src+K;
			f[3][1] = v1;
			
			//分支节点4
			if(state==1)v1=0;
			else v1=Math.abs(state-1)+K;
			f[4][0] = v1;
			
			if(state==0) v1=0;
			else v1= Math.abs(state)+K;
			f[4][1] = v1;
			
			//分支节点5
			if(p==0||tag==9999)v1=0;
			else v1=Math.min(Math.abs(p-0)+K,Math.abs(tag-9999)+K);
			f[5][0] = v1;
			
			if(p!=0&&tag!=9999)v1=0;
			else v1= 2*K;
			f[5][1] = v1;
		}
		if(func_num == 4)
		{
			int Direction = (int)x[0];
			String map1 = String.valueOf((char)x[1])+String.valueOf((char)x[2])+String.valueOf((char)x[3]);
			String map2 = String.valueOf((char)x[4])+String.valueOf((char)x[5])+String.valueOf((char)x[6]);
			
			double v1=0;
			
			//分支节点1
			if(Direction==1) v1=0;
			else v1 = Math.abs(Direction-1)+K;
			f[0][0] = v1;
			
			if(Direction!=1) v1=0;
			else v1 = K;
			f[0][1] = v1;
			
			//分支节点2
			if(map1.equals("001")||map1.equals("002")||map1.equals("003"))
				v1=0;
			else
//				v1 = Math.abs((char)x[1]-'0') + Math.abs((char)x[2]-'0') + Math.min(Math.min((char)x[3]-'1', (char)x[3]-'2'), (char)x[3]-'3');
			{
				v1 = Math.min(Math.abs((char)x[3]-'1')+K, Math.abs((char)x[3]-'2')+K);
				v1 = Math.min(v1, Math.abs((char)x[3]-'3')+K);
				v1 = Math.abs((char)x[1]-'0') + Math.abs((char)x[2]-'0') + v1;
			}
			f[1][0] = v1;
			
			if(!map1.equals("001")&&!map1.equals("002")&&!map1.equals("003"))
				v1 = 0;
			else
				v1 = 3*K;
			f[1][1] = v1;
			
			//分支节点3
			if(map2.equals("004")||map2.equals("005")||map2.equals("006"))
				v1=0;
			else
//				v1 = Math.abs((char)x[4]-'0') + Math.abs((char)x[5]-'0') + Math.min(Math.min((char)x[6]-'4', (char)x[6]-'5'), (char)x[6]-'6');
			{
				v1 = Math.min(Math.abs((char)x[6]-'4')+K, Math.abs((char)x[6]-'5')+K);
				v1 = Math.min(v1, Math.abs((char)x[6]-'6')+K);
				v1 = Math.abs((char)x[4]-'0') + Math.abs((char)x[5]-'0') + v1;
			}
			f[2][0] = v1;
			
			if(!map2.equals("004")&&!map2.equals("005")&&!map2.equals("006"))
				v1 = 0;
			else
				v1 = 3*K;
			f[2][1] = v1;
		}
		if(func_num == 5)
		{
			String cl = String.valueOf((char)x[1])+String.valueOf((char)x[2])+String.valueOf((char)x[3]);
			int p1 = (int)x[0], p2 = (int)x[R-1];
			
			double v1=0;
			//分支节点1
			if(p1==0) v1=0;
			else v1=Math.abs(p1)+K;
			f[0][0] = v1;
			
			if(p1==1) v1=0;
			else v1 = Math.abs(p1-1)+K;
			f[0][1] = v1;
			
			//分支节点2
			if(cl.equals("001")||cl.equals("002")||cl.equals("003"))
				v1=0;
			else
//				v1=Math.abs((char)x[1]-'0') + Math.abs((char)x[2]-'0') + Math.min(Math.min((char)x[3]-'1', (char)x[3]-'2'), (char)x[3]-'3');
			{
				v1 = Math.min(Math.abs((char)x[3]-'1')+K, Math.abs((char)x[3]-'2')+K);
				v1 = Math.min(v1, Math.abs((char)x[3]-'3')+K);
				v1 = Math.abs((char)x[1]-'0') + Math.abs((char)x[2]-'0') + v1;
			}
			f[1][0] = v1;
			
			if(!cl.equals("001")&&!cl.equals("002")&&!cl.equals("003"))
				v1 = 0;
			else
				v1 = 3*K;
			f[1][1] = v1;
			
			//分支节点3
			if(p2==0) v1=0;
			else v1=Math.abs(p2)+K;
			f[2][0] = v1;
			
			if(p2==1) v1=0;
			else v1 = Math.abs(p2-1)+K;
			f[2][1] = v1;
		}
		if(func_num == 6)
		{
			String edge = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2]); 
			String pair = String.valueOf((char)x[3])+String.valueOf((char)x[4])+String.valueOf((char)x[5]); 
			
			double v1=0;
			//分支节点1
			if(edge.equals("mod")) v1=0;
			else v1 = Math.abs((char)x[0]-'m')+Math.abs((char)x[1]-'o')+Math.abs((char)x[2]-'d')+3*K;
			f[0][0] = v1;
			
			if(!edge.equals("mod")) v1=0;
			else v1 = K;
			f[0][1] = v1;
			
			//分支节点2
			if(pair.equals("001")||pair.equals("002")||pair.equals("003"))
				v1 = 0;
			else
//				v1 = Math.abs((char)x[3]-'0') + Math.abs((char)x[4]-'0') + Math.min(Math.min((char)x[5]-'1', (char)x[5]-'2'), (char)x[5]-'3');
			{
				v1 = Math.min(Math.abs((char)x[5]-'1')+K, Math.abs((char)x[5]-'2')+K);
				v1 = Math.min(v1, Math.abs((char)x[5]-'3')+K);
				v1 = Math.abs((char)x[3]-'0') + Math.abs((char)x[4]-'0') + v1;
			}
			f[1][0] = v1;
			
			if(!pair.equals("001")&&!pair.equals("002")||!pair.equals("003"))
				v1 = 0;
			else
				v1 = 3*K;
			f[1][1] = v1;
			//分支节点3
			if((int)x[6]==0) v1 = 0;
			else v1 = Math.abs((int)x[6])+K;
			f[2][0] = v1;
			
			if((int)x[6]==1) v1 =0;
			else v1 = Math.abs((int)x[6]-1)+K;
			f[2][1] = v1;
			//分支节点4
			if((int)x[7]==2) v1=0;
			else v1 = Math.abs((int)x[7]-2)+K;
			f[3][0] = v1;
			
			if((int)x[7]!=2) v1=0;
			else v1 = K;
			f[3][1] = K;
		}

		if(func_num==7){
			double r,g,b,max,min,h;
			r = x[0]; g = x[1]; b = x[2];
			max = Math.max(Math.max(r, g ), b);
			min = Math.min(Math.min(r, g ), b);
			h = max-min>0.0001&&r==max?(g-b)/(max-min):0;
			double v1;
			//分支节点1
			if(max-min>0.0001) v1=0;
			else v1 = Math.abs(max-min-0.0001)+K;
			f[0][0] = v1;
			
			if(max-min<=0.0001) v1=0;
			else v1 = Math.abs(max-min-0.0001)+K;
			f[0][1] = v1;
			//分支节点2
			if(r==max) v1=0;
			else v1 = Math.abs(r-max)+K;
			f[1][0] = v1;
			
			if(r!=max) v1=0;
			else v1 = K;
			f[1][1] = v1;
			//分支节点3
			if(h<0.0) v1=0;
			else v1 = h+K;
			f[2][0] = v1;
			
			if(h>=0.0) v1=0;
			else v1 = -h+K;
			f[2][1] = v1;
			//分支节点4
			if(g==max) v1=0;
			else v1=Math.abs(g-max)+K;
			f[3][0] = v1;
			
			if(g!=max) v1=0;
			else v1=K;
			f[3][0] = v1;
			//分支节点5
			if(b==max) v1=0;
			else v1= Math.abs(b-max)+K;
			f[4][0] = v1;
			
			if(b!=max) v1=0;
			else v1=K;
			f[4][1] = v1;
		}
		if(func_num==8){
			double h,s,v,hue;
			int i;
			h=x[0];s=x[1];v=x[2];
			hue = s==0?0:h;
			double v1;
			//分支节点1
			if((int)s==0) v1=0;
			else v1=Math.abs((int)s)+K;
			f[0][0] = v1;
			
			if((int)s!=0) v1=0;
			else v1 = K;
			f[0][1] = v1;
			//分支节点2
			if((int)hue==1) v1=0;
			else v1=Math.abs((int)hue-1)+K;
			f[1][0] = v1;
			
			if((int)hue!=1) v1=0;
			else v1=K;
			f[1][1] = v1;
			//分支节点3
			i = (int) ( 6*((int)hue==1?0:hue));
			if(i==0)	f[2][0]=0;
			else	f[2][0] = Math.abs(i-0)+K;
			if(i==1)	f[2][1]=0;
			else	f[2][1] = Math.abs(i-1)+K;
			if(i==2)	f[2][2]=0;
			else	f[2][2] = Math.abs(i-2)+K;
			if(i==3)	f[2][3]=0;
			else	f[2][3] = Math.abs(i-3)+K;
			if(i==4)	f[2][4]=0;
			else	f[2][4] = Math.abs(i-4)+K;
			if(i==5)	f[2][5]=0;
			else	f[2][5] = Math.abs(i-5)+K;
			if(i==6)	f[2][6]=0;
			else	f[2][6] = Math.abs(i-6)+K;
		}
		if(func_num==9){
			String ISBN = null; int k= (int) x[11];
			
			for(int i=0;i<R-1;i++)
				ISBN += (char)x[i];
			double v1;
			
			//分支节点1
			if(ISBN.charAt(k)==' '||ISBN.charAt(k)=='-') {f[0][0]=0; f[0][1]=Math.min(Math.abs(ISBN.charAt(k)-' '), Math.abs(ISBN.charAt(k)-'-'))+K;}
			else if((ISBN.charAt(k)-'0'>=0&&ISBN.charAt(k)-'0'<=9)||ISBN.charAt(k)=='X'||ISBN.charAt(k)=='x'){
				f[0][0] = 0;
				f[0][1] = Math.min(Math.abs(ISBN.charAt(k)-'X'), Math.abs(ISBN.charAt(k)-'x'));
				f[0][1] = Math.min(f[0][1], Math.abs(ISBN.charAt(k)-'0'));
				f[0][1] = Math.min(f[0][1], Math.abs(ISBN.charAt(k)-'1'));
				f[0][1] = Math.min(f[0][1], Math.abs(ISBN.charAt(k)-'2'));
				f[0][1] = Math.min(f[0][1], Math.abs(ISBN.charAt(k)-'3'));
				f[0][1] = Math.min(f[0][1], Math.abs(ISBN.charAt(k)-'4'));
				f[0][1] = Math.min(f[0][1], Math.abs(ISBN.charAt(k)-'5'));
				f[0][1] = Math.min(f[0][1], Math.abs(ISBN.charAt(k)-'6'));
				f[0][1] = Math.min(f[0][1], Math.abs(ISBN.charAt(k)-'7'));
				f[0][1] = Math.min(f[0][1], Math.abs(ISBN.charAt(k)-'8'));
				f[0][1] = Math.min(f[0][1], Math.abs(ISBN.charAt(k)-'9'));
				f[0][1] += K;
			}else {f[0][0] = K; f[0][1] = 0;}
			//分支节点2
			if(k<10) v1=0;
			else v1 = Math.abs(k-10)+K;
			f[1][0] = v1;
			
			if(k>=10) v1=0;
			else v1 = Math.abs(k-10)+K;
			f[1][1] = v1;
			//分支节点3
			if(k==10) v1=0;
			else v1=Math.abs(k-10)+K;
			f[2][0] = v1;
			
			if(k!=10) v1=0;
			else v1 = K;
			f[2][1] = v1;
			//分支节点4
			int temp = (ISBN.charAt(10)-'X'==0||ISBN.charAt(10)-'x'==0)?10:0;
			if(temp!=10)
				temp = (ISBN.charAt(10)-'0'>=0&&ISBN.charAt(10)-'0'<=9)?ISBN.charAt(10)-'0':0;
			if(checksum(ISBN)%11!=temp) v1=0;
			else v1 = K;
			f[3][0] = v1;
			
			if(checksum(ISBN)%11==temp) v1=0;
			else v1 = Math.abs(checksum(ISBN)%11-temp)+K;
			f[3][1] = v1;
			//分支节点5
			if(k>0) v1=0;
			else v1 = Math.abs(k)+K;
			f[4][0] = v1;
			
			if(k<=0) v1=0;
			else v1 = Math.abs(k)+K;
			f[4][1] = v1;
		}
		if(func_num==10){
			String ISSN = null; int k= (int) x[9];
			if(k>=ub[R-1]) k = k-1;
			for(int i=0;i<9;i++)
				ISSN += (char)x[i];
			double v1;
			
			//分支节点1
			if(ISSN.charAt(k)==' '||ISSN.charAt(k)=='-') {f[0][0]=0; f[0][1]=Math.min(Math.abs(ISSN.charAt(k)-' '), Math.abs(ISSN.charAt(k)-'-'))+K;}
			else if((ISSN.charAt(k)-'0'>=0&&ISSN.charAt(k)-'0'<=9)||ISSN.charAt(k)=='X'||ISSN.charAt(k)=='x'){
				f[0][0] = 0;
				f[0][1] = Math.min(Math.abs(ISSN.charAt(k)-'X'), Math.abs(ISSN.charAt(k)-'x'));
				f[0][1] = Math.min(f[0][1], Math.abs(ISSN.charAt(k)-'0'));
				f[0][1] = Math.min(f[0][1], Math.abs(ISSN.charAt(k)-'1'));
				f[0][1] = Math.min(f[0][1], Math.abs(ISSN.charAt(k)-'2'));
				f[0][1] = Math.min(f[0][1], Math.abs(ISSN.charAt(k)-'3'));
				f[0][1] = Math.min(f[0][1], Math.abs(ISSN.charAt(k)-'4'));
				f[0][1] = Math.min(f[0][1], Math.abs(ISSN.charAt(k)-'5'));
				f[0][1] = Math.min(f[0][1], Math.abs(ISSN.charAt(k)-'6'));
				f[0][1] = Math.min(f[0][1], Math.abs(ISSN.charAt(k)-'7'));
				f[0][1] = Math.min(f[0][1], Math.abs(ISSN.charAt(k)-'8'));
				f[0][1] = Math.min(f[0][1], Math.abs(ISSN.charAt(k)-'9'));
				f[0][1] += K;
			}else {f[0][0] = K; f[0][1] = 0;}
			//分支节点2
			if(k<8) v1=0;
			else v1 = Math.abs(k-8)+K;
			f[1][0] = v1;
			
			if(k>=8) v1=0;
			else v1 = Math.abs(k-8)+K;
			f[1][1] = v1;
			//分支节点3
			if(k==8) v1=0;
			else v1=Math.abs(k-8)+K;
			f[2][0] = v1;
			
			if(k!=8) v1=0;
			else v1 = K;
			f[2][1] = v1;
			//分支节点4
			int temp = (ISSN.charAt(8)-'X'==0||ISSN.charAt(8)-'x'==0)?10:0;
			if(temp!=10)
				temp = (ISSN.charAt(8)-'0'>=0&&ISSN.charAt(8)-'0'<=9)?ISSN.charAt(8)-'0':0;
			if(checksum(ISSN)%11!=temp) v1=0;
			else v1 = K;
			f[3][0] = v1;
			
			if(checksum(ISSN)%11==temp) v1=0;
			else v1 = Math.abs(checksum(ISSN)%11-temp)+K;
			f[3][1] = v1;
			//分支节点5
			if(k>0) v1=0;
			else v1 = Math.abs(k)+K;
			f[4][0] = v1;
			
			if(k<=0) v1=0;
			else v1 = Math.abs(k)+K;
			f[4][1] = v1;
		}
		if(func_num==11){
			double h =  x[0];double w = x[1];double b = x[2];
			int i;
			h = 6.0* h/360.0;
			double temp = Math.floor(h);
			i = (int) temp;
			double v1;
			
			//分支节点1
			if(temp==-1.0) v1=0;
			else v1 = Math.abs(temp+1.0)+ K;
			f[0][0] = v1;
			
			if(temp!=-1.0) v1=0;
			else v1 = K;
			f[0][1] = v1;
			
			//分支节点2
			if(temp==1.0) v1=0;
			else v1=Math.abs(temp-1.0)+K;
			f[1][0] = v1;
			
			if(temp!=1.0) v1=0;
			else v1=K;
			f[1][1] = v1;
			//分支节点3
			if(i==0)	f[2][0]=0;
			else	f[2][0] = Math.abs(i)+K;
			if(i==1)	f[2][1]=0;
			else	f[2][1] = Math.abs(i-1)+K;
			if(i==2)	f[2][2]=0;
			else	f[2][2] = Math.abs(i-2)+K;
			if(i==3)	f[2][3]=0;
			else	f[2][3] = Math.abs(i-3)+K;
			if(i==4)	f[2][4]=0;
			else	f[2][4] = Math.abs(i-4)+K;
			if(i==5)	f[2][5]=0;
			else	f[2][5] = Math.abs(i-5)+K;
			if(i==6)	f[2][6]=0;
			else	f[2][6] = Math.abs(i-6)+K;
		}
		if(func_num==12){
			int nSample = (int)x[0], nPixel = (int) x[1];
			double v1;
			
			//分支节点1
			if(nSample==1&&nPixel==1) v1=0;
			else v1 = Math.abs(nSample-1)+Math.abs(nPixel-1)+2*K;
			f[0][1] = v1;
			
			if(nSample!=1||nPixel!=1) v1=0;
			else v1 = K;
			f[0][1] = v1;
			//分支节点2
			if(nSample==1&&nPixel==2) v1=0;
			else v1 = Math.abs(nSample-1)+Math.abs(nPixel-2)+2*K;
			f[1][0] = v1;
			
			if(nSample!=1||nPixel!=2) v1=0;
			else v1 = K;
			f[1][1] = v1;
			//分支节点3
			if(nSample==1&&nPixel==4) v1=0;
			else v1 = Math.abs(nSample-1)+Math.abs(nPixel-4)+2*K;
			f[2][0] = v1;
			
			if(nSample!=1||nPixel!=4) v1=0;
			else v1 = K;
			f[2][1] = v1;
			//分支节点4
			if(nSample==2&&nPixel==2) v1=0;
			else v1 = Math.abs(nSample-2)+Math.abs(nPixel-2)+2*K;
			f[3][0] = v1;
			
			if(nSample!=2||nPixel!=2) v1=0;
			else v1 = K;
			f[3][1] = v1;
			//分支节点5
			if(nSample==2&&nPixel==32) v1=0;
			else v1 = Math.abs(nSample-2)+Math.abs(nPixel-32)+2*K;
			f[4][0] = v1;
			
			if(nSample!=2||nPixel!=32) v1=0;
			else v1 = K;
			f[4][1] = v1;
			//分支节点6
			if(nSample==3&&nPixel==4) v1=0;
			else v1 = Math.abs(nSample-3)+Math.abs(nPixel-4)+2*K;
			f[5][1] = v1;
			
			if(nSample!=3||nPixel!=4) v1=0;
			else v1 = K;
			f[5][1] = v1;
			//分支节点7
			if(nSample==3&&nPixel==8) v1=0;
			else v1 = Math.abs(nSample-3)+Math.abs(nPixel-8)+2*K;
			f[6][0] = v1;
			
			if(nSample!=3||nPixel!=8) v1=0;
			else v1 = K;
			f[6][1] = v1;
		}
		
		if(func_num == 13)
		{
			String option = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2]); 
			String extraOption = String.valueOf((char)x[3])+String.valueOf((char)x[4])+String.valueOf((char)x[5]); 
			int type = x[6];
			
			double v1=0;
			//分支节点1
			if(option.equals("   ")) v1 = 0;
			else v1 = Math.abs((char)x[0]-' ')+Math.abs((char)x[1]-' ')+Math.abs((char)x[2]-' ')+3*K;
			f[0][0] = v1;
			
			if(!option.equals("   ")) v1 =0;
			else v1 = K;
			f[0][1] = v1;
			
			//分支节点2
			if(!extraOption.equals("   ")) v1=0;
			else v1 = K;
			f[1][0] = v1;
			
			if(extraOption.equals("   ")) v1=0;
			else v1 = Math.abs((char)x[3]-' ')+Math.abs((char)x[4]-' ')+Math.abs((char)x[5]-' ')+3*K;
			f[1][1] = v1;
			
			//分支节点3
			if(x[5]==',') v1=0;
			else v1 = Math.abs((char)x[5] - ',') + K;
			f[2][0] = v1;
			
			if(x[5]!=',') v1 = 0;
			else v1 = K;
			f[2][1] = v1;
			
			//分支节点4
			if(type==0)	f[3][0] = 0;
			else	f[3][0] = Math.abs(type-0)+K;
			if(type==1)	f[3][1] = 0;
			else	f[3][1] = Math.abs(type-1)+K;
			if(type==2)	f[3][2] = 0;
			else	f[3][2] = Math.abs(type-2)+K;
			if(type==3)	f[3][3] = 0;
			else	f[3][3] = Math.abs(type-3)+K;
			if(type==4)	f[3][4] = 0;
			else	f[3][4] = Math.abs(type-4)+K;
			if(type==5)	f[3][5] = 0;
			else	f[3][5] = Math.abs(type-5)+K;
			if(type==6)	f[3][6] = 0;
			else	f[3][6] = Math.abs(type-6)+K;
		}
		if(func_num == 14){
			String xmlTags = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2]); 
			String sentence = String.valueOf((char)x[3])+String.valueOf((char)x[4])+String.valueOf((char)x[5]); 
			
			double v1=0;
			//分支节点1
			if(xmlTags.equals("   ")) v1 = 0;
			else v1 = Math.abs((char)x[0]-' ')+Math.abs((char)x[1]-' ')+Math.abs((char)x[2]-' ')+3*K;
			f[0][0] = v1;
			
			if(!xmlTags.equals("   ")) v1 =0;
			else v1 = K;
			f[0][1] = v1;
			
			//分支节点2
			if(sentence.equals("   ")) v1 = 0;
			else v1 = Math.abs((char)x[3]-' ')+Math.abs((char)x[4]-' ')+Math.abs((char)x[5]-' ')+3*K;
			f[1][0] = v1;
			
			if(!sentence.equals("   ")) v1 =0;
			else v1 = K;
			f[1][1] = v1;			
		}
		if(func_num == 15){
			boolean nlSplitting,whitespaceTokenization;
			if(x[0]==0) nlSplitting = true;
			else nlSplitting = false;
			if(x[1]==0) whitespaceTokenization = true;
			else whitespaceTokenization = false;
			String line = String.valueOf((char)x[2])+String.valueOf((char)x[3]); 
			String isOneSentence = String.valueOf((char)x[4])+String.valueOf((char)x[5])+String.valueOf((char)x[6])+String.valueOf((char)x[7]); 
			char token,bound1,bound2;
			token = (char)x[8]; bound1 = (char)x[9]; bound2 = (char)x[10];
			
			double v1=0;
			//分支节点1
			if(nlSplitting) v1 = 0;
			else v1 = K;
			f[0][0] = v1;
			
			if(!nlSplitting) v1 = 0;
			else v1 = K;
			f[0][1] = v1;
			
			//分支节点2
			if(whitespaceTokenization) v1 = 0;
			else v1 = K;
			f[1][0] = v1;
			
			if(!whitespaceTokenization) v1 = 0;
			else v1 = K;
			f[1][1] = v1;
			
			//分支节点3
			if(line.equals("/n")) v1 = 0;
			else v1 = Math.abs((char)x[2]-'/')+Math.abs((char)x[3]-'n')+2*K;
			f[2][0] = v1;
			
			if(!line.equals("\n")) v1 = 0;
			else v1 = K;
			f[2][1] = v1;
			
			//分支节点4
			if(isOneSentence.equals("true")) v1 = 0;
			else v1 = Math.abs((char)x[4]-'t')+Math.abs((char)x[5]-'r')+Math.abs((char)x[6]-'u')+Math.abs((char)x[7]-'e')+4*K;
			f[3][0] = v1;
			
			if(!isOneSentence.equals("true")) v1 = 0;
			else v1 = K;
			f[3][1] = v1;
			
			//分支节点5
			if(token==' ') v1 = 0;
			else v1 = Math.abs(token-' ')+K;
			f[4][0] = v1;
			
			if(token!=' ') v1 = 0;
			else v1 = K;
			f[4][1] = v1;
			
			//分支节点6
			if(bound1==' ') v1 = 0;
			else v1 = Math.abs(bound1-' ')+K;
			f[5][0] = v1;
			
			if(bound1!=' ') v1 = 0;
			else v1 = K;
			f[5][1] = v1;
			
			//分支节点7
			if(bound2==' ') v1 = 0;
			else v1 = Math.abs(bound2-' ')+K;
			f[6][0] = v1;
			
			if(bound2!=' ') v1 = 0;
			else v1 = K;
			f[6][1] = v1;
		}
		if(func_num == 16){
			String annotation = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2]); 
			int nThreads = x[3];
			
			double v1;
			//分支节点1
			if(annotation.equals("001")||annotation.equals("002")||annotation.equals("003")) v1 = 0;
			else {
				v1 = Math.min(Math.abs((char)x[2]-'1')+K, Math.abs((char)x[2]-'2')+K);
				v1 = Math.min(v1, Math.abs((char)x[2]-'3')+K);
				v1 = Math.abs((char)x[0]-'0') + Math.abs((char)x[1]-'0') + v1;
			}
			f[0][0] = v1;
			
			if(!annotation.equals("001")&&!annotation.equals("002")||!annotation.equals("003"))
				v1 = 0;
			else
				v1 = 3*K;
			f[0][1] = v1;
			//分支节点2
			if(nThreads == 1) v1 = 0;
			else v1 = Math.abs(nThreads-1)+K;
			f[1][0] = v1;
			
			if(nThreads != 1) v1=0;
			else v1 = K;
			f[1][1] = v1;
		}
		if(func_num == 17){
			String nerLanguage = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2])+
					String.valueOf((char)x[3])+String.valueOf((char)x[4])+String.valueOf((char)x[5])+String.valueOf((char)x[6]); 
			boolean augment;
			if(x[7]==0) augment = true;
			else augment = false;
			
			double v1;
			//分支节点1
			if(nerLanguage.equals("CHINESE")) v1 = 0;
			else v1 = Math.abs((char)x[0]-'C')+Math.abs((char)x[1]-'H')+Math.abs((char)x[2]-'I')+Math.abs((char)x[3]-'N')+
					Math.abs((char)x[4]-'E')+Math.abs((char)x[5]-'S')+Math.abs((char)x[6]-'E')+7*K;
			f[0][0] = v1;
			
			if(!nerLanguage.equals("CHINESE")) v1 = 0;
			else v1 = 7*K;
			f[0][1] = v1;
			
			//分支节点2
			if(augment) v1=0;
			else v1 = K;
			f[1][0] = v1;
			
			if(!augment) v1=0;
			else v1 = K;
			f[1][1] = v1;
		}
		if(func_num == 18){
			String trueCase = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2])+String.valueOf((char)x[3])+String.valueOf((char)x[4]);
			boolean overwriteText;
			if(x[5]==0) overwriteText = true;
			else overwriteText = false;
			
//			double v1;
			//分支节点1
			if(trueCase.equals("UPPER")) 
				f[0][0] = 0;
			else 
				f[0][0] = Math.abs((char)x[0]-'U')+Math.abs((char)x[1]-'P')+Math.abs((char)x[2]-'P')+Math.abs((char)x[3]-'E')+Math.abs((char)x[4]-'R')+5*K;

			if(!trueCase.equals("UPPER")) 
				f[0][1] = 0;
			else
				f[0][1] = 5*K;
			
			//分支节点2
			if(trueCase.equals("LOWER")) 
				f[1][0] = 0;
			else 
				f[1][0] = Math.abs((char)x[0]-'L')+Math.abs((char)x[1]-'O')+Math.abs((char)x[2]-'W')+Math.abs((char)x[3]-'E')+Math.abs((char)x[4]-'R')+5*K;

			if(!trueCase.equals("LOWER")) 
				f[1][1] = 0;
			else
				f[1][1] = 5*K;
			
			//分支节点3
			if(trueCase.equals("INIT_"))
				f[2][0] = 0;
			else
				f[2][0] = Math.abs((char)x[0]-'I')+Math.abs((char)x[1]-'N')+Math.abs((char)x[2]-'I')+Math.abs((char)x[3]-'T')+Math.abs((char)x[4]-'_')+5*K;
			
			if(!trueCase.equals("INIT_")) 
				f[2][1] = 0;
			else
				f[2][1] = 5*K;
			
			//分支节点4
			if(trueCase.equals("O    "))
				f[3][0] = 0;
			else
				f[3][0] = Math.abs((char)x[0]-'O')+Math.abs((char)x[1]-' ')+Math.abs((char)x[2]-' ')+Math.abs((char)x[3]-' ')+Math.abs((char)x[4]-' ')+5*K;
			
			if(!trueCase.equals("O   ")) 
				f[3][1] = 0;
			else
				f[3][1] = 5*K;
			
			//分支节点5
			if(overwriteText) {f[4][0] = 0; f[4][1] = K;}
			else {f[4][0] = K; f[4][1] = 0;}
		}
		if(func_num == 19)   //Triangle
		{
			int a = x[0] ;
			int b = x[1] ;
			int c = x[2] ;		
					
    		double v1,v2,v3,v4,v5,v6,v7 ;	   		   	   		
	   		
    		if(a<(b+c)) v1 = 0;        //测试用例执行第一个节点的Yes分支时的成本值f[0]；
	   		else v1 = a-(b+c)+K ;
	   		if(b<(a+c)) v2 = 0 ;
	   		else v2 = b-(a+c) + K ;
	   		if(c<(a+b)) v3 = 0;
	   		else v3 = c-(a+b) + K ;	  		
	   		f[0][0] = v1 + v2 + v3 ;  
	   		
	   		if(a==b) v1 = 0 ;      //测试用例执行第二个节点的Yes分支时的成本值f[1]；
	   		else v1 = Math.abs(a-b)+K ;
	   		if(a!=c) v2 = 0 ;
	   		else v2 = K ;
	   		if(a==c) v3 = 0 ;
	   		else v3 = Math.abs(a-c)+K ;
	   		if(a!=b) v4 = 0 ;
	   		else v4 = K ;
	   		if(b==c) v5 = 0 ;
	   		else v5 = Math.abs(b-c)+K ;
	   		if(b!=a) v6 = 0 ;
	   		else v6 = K ;
	   		v7 = Math.min(v1+v2 , v3+v4);
	   		f[1][0] = Math.min(v7 , v5+v6);
	   		
	   		if(a==b) v1 = 0 ;     //测试用例执行第三个节点的Yes分支时的成本值f[2]；
	   		else v1 = Math.abs(a-b)+K ;
	   		if(a==c) v2 = 0;
	   		else v2 = Math.abs(a-c)+K ;
	   		f[2][0] = v1 + v2 ;
	   		
	   		if(a!=b) v1 = 0 ;   //测试用例执行第四个节点的Yes分支时的成本值f[3]；
	   		else v1 = K ;
	   		if(a!=c) v2 = 0 ;
	   		else v2 = K ;
	   		if(b!=c) v3 = 0 ;
	   		else v3 = K ;
	   		f[3][0] = v1 + v2 + v3 ;

	   		if(a>=(b+c)) v1 = 0 ;   //测试用例执行第一个节点的No分支时的成本值F[0]；
	   		else v1 = (b+c)-a+K ;
	   		if(b>=(a+c)) v2 = 0;
	   		else v2 = (a+c)-b+K ;
	   		if(c>=(a+b)) v3 = 0 ;
	   		else v3 = (a+b)-c+K ;
	   		v4 = Math.min(v1, v2);
	   		f[0][1] = Math.min(v4, v3);
   		
	   		if(a!=b) v1 = 0 ;       //测试用例执行第二个节点的No分支时的成本值F[1]；
	   		else v1 = K ;
	   		if(a==c) v2 = 0 ;
	   		else v2 = Math.abs(a-c)+K ;
	   		if(a!=c) v3 = 0 ;
	   		else v3 = K ;
	   		if(a==b) v4 = 0 ;
	   		else v4 = Math.abs(a-b)+K ;
	   		if(b!=c) v5 = 0 ;
	   		else v5 = K ;
	   		if(b==a) v6 = 0 ;
	   		else v6 = Math.abs(b-a)+K ;
	   		f[1][1] = Math.min(v1, v2) + Math.min(v3, v4) + Math.min(v5, v6) ;
  		
	   		if(a!=b) v1 = 0 ;     //测试用例执行第三个节点的No分支时的成本值F[2]；
	   		else v1 = K ;
	   		if(a!=c) v2 = 0 ;
	   		else v2 = K ;
	   		f[2][1] = Math.min(v1, v2);
   		
	   		if(a==b) v1 = 0 ;    //测试用例执行第四个节点的No分支时的成本值F[3]；
	   		else v1 = Math.abs(a-b)+K ;
	   		if(a==c) v2 = 0 ;
	   		else v2 = Math.abs(a-c)+K ;
	   		if(b==c) v3 = 0 ;
	   		else v3 = Math.abs(b-c)+K ;
	   		v4 = Math.min(v1, v2);
	   		f[3][1] = Math.min(v4, v3);		
		}	  
	    if(func_num == 20)    //Factorial
	   	{
	   		int a = x[0];
	   		double v1 ;
	   		
	   		if(a==1) v1 = 0 ;
	   		else v1 = Math.abs(a-1)+K ;
	   		f[0][0] = v1 ;
	   		  
	   		if(a!=1) v1 = 0 ;
	   		else v1 = K ;
	   		f[0][1] = v1 ;
	   	 }
	 	  if(func_num == 21)   //sorting
	   	  {
	   		  boolean d =false ;
	   		  double v1=0,v2=0 ;
		   	  int i1,j1 ;
		   	  int[] a = new int[R];
		   	  for(i1=0;i1<R;i1++)
		   	       a[i1] = x[i1];

		   	  for(j1=0;j1<=R-1;j1++) 
		   	  {
		   		  for (i1=0;i1<R-1-j1;i1++)
		   		  {
		   			  d =(a[i1]>a[i1+1]) ;
		   			  if(a[i1]>a[i1+1]){ v1 = 0 ; break ;}
		   			else v1 = a[i1+1]-a[i1]+K ;
		   			  v2 = v2 + v1 ;
		   		  }
		   		  if(d) break;	   		  
		   	  }
		   	  if(v1==0) f[0][0] = v1 ;
		   	  else f[0][0] = v2 ;
		   	 
		   	  v2 = 0 ;
			  for(j1=0;j1<=R-1;j1++) 
		   	  {
		   		  for (i1=0;i1<R-1-j1;i1++)
		   		  {
		   			  if(a[i1]<=a[i1+1]) v1 = 0 ; 
		   			else v1 = a[i1]-a[i1+1]+K ;
		   			  v2 = v2 + v1 ;
		   		  }
		   	  }
		   	  f[0][1] = v2 ;		   		
	   	  }
	 	  if(func_num == 22)    //GCD
	   	  {
	   		 int m = x[0] ;
	   		 int n = x[1] ;
	   		 double v1 ;
	   		 	   		
	   		 if(m<n)v1 = 0 ;
	   		 else v1 = m-n+K ;
	   		 f[0][0] = v1 ;
	   		
	   		 if(m>=n)v1 = 0 ;
	   		 else v1 = n-m+K ;
	   		 f[0][1] = v1 ;
	   		
	   	     int r;
	   	     r = m % n;
	   		 m = n;
	   		 n = r;
	   		
	   		 if(r!=0) v1 = 0 ;
	   		else v1 = K ;
	   		 f[1][0] = v1 ;
	   		
	   		 if(r==0)v1 = 0 ;
	   		else v1 = Math.abs(r)+K ;
	   		 f[1][1] = v1 ;
	   	  }
	 	  
	 	 if(func_num == 23)    //Middle
	     {
    	     int a = x[0] ;
	   		 int b = x[1] ;
	   		 int c = x[2] ;
	   		 double v1,v2,v3,v4 ;		   		
	   		
	   		 if(a<b) v1 = 0 ;
	   		 else v1 = a-b+K ;
	   		 if(b<c) v2 = 0 ;
	   		 else v2 = b-c+K ;
	   		 if(c<b) v3 = 0 ;
	   		 else v3 = c-b+K;
	   		 if(b<a) v4 = 0 ;
	   		 else v4 = b-a+K ;
	   		 f[0][0] = Math.min(v1+v2, v3+v4);
	   		
	   		 if(a>=b) v1 = 0 ;
	   		 else v1 = b-a+K ;
	   		 if(b>=c) v2 = 0 ;
	   		 else v2 = c-b+K ;
	   		 if(c>=b) v3 = 0 ;
	   		 else v3 = b-c+K;
	   		 if(b>=a) v4 = 0 ;
	   		 else v4 = a-b+K ;
	   		 f[0][1] = Math.min(v1, v2) + Math.min(v3, v4);
	   		
	   		 if(a<c) v1 = 0 ;
	   		 else v1 = a-c+K ;
	   		 if(c<b) v2 = 0 ;
	   		 else v2 = c-b+K ;
	   		 if(b<c) v3 = 0 ;
	   		 else v3 = b-c+K;
	   	 	 if(c<a) v4 = 0 ;
	   		 else v4 = c-a+K ;
	   		 f[1][0] = Math.min(v1+v2, v3+v4);
	   		
	   		 if(a>=c) v1 = 0 ;
	   		 else v1 = c-a+K ;
	   		 if(c>=b) v2 = 0 ;
	   		 else v2 = b-c+K ;
	   		 if(b>=c) v3 = 0 ;
	   		 else v3 = c-b+K;
	   		 if(c>=a) v4 = 0 ;
	   		 else v4 = a-c+K ;
	   		 f[1][1] = Math.min(v1, v2) + Math.min(v3, v4);
	   		
	   		 if(b<a) v1 = 0 ;
	   		 else v1 = b-a+K ;
	   		 if(a<c) v2 = 0 ;
	   		 else v2 = a-c+K ;
	   		 if(c<a) v3 = 0 ;
	   		 else v3 = c-a+K;
	   		 if(a<b) v4 = 0 ;
	   		 else v4 = a-b+K ;
	   		 f[2][0] = Math.min(v1+v2, v3+v4);
	   		
	   		 if(b>=a) v1 = 0 ;
	   		 else v1 = a-b+K ;
	   		 if(a>=c) v2 = 0 ;
	   		 else v2 = c-a+K ;
	   		 if(c>=a) v3 = 0 ;
	   		 else v3 = a-c+K;
	   		 if(a>=b) v4 = 0 ;
	   		 else v4 = b-a+K ;
	   		 f[2][1] = Math.min(v1, v2) + Math.min(v3, v4);
	     }
	 	if(func_num == 24)   //Tomorrow
      	{
      		int Day = x[0] ;
      	    int Year = x[1] ;
      	    int Month = x[2] ;
      	    int Date = x[3] ;
		   		
	   		double v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15 ;	   		
	   		   	    
	   		if(Day == 7) v1 = 0 ;
	   		else v1 = Math.abs(Day-7)+K ;
	   		f[0][0] = v1 ;
	   		
	   		if(Day != 7) v1 = 0 ;
	   		else v1 = K ;
	   		f[0][1] = v1 ;

	   	    if(Month == 12) v1 = 0 ;
	   	    else v1 = Math.abs(Month-12)+K ;
	   	    if(Date == 31) v2 = 0 ;
	   	    else v2 = Math.abs(Date-31)+K;
	   	    f[1][0] = v1 + v2 ;
	   	     
	   	    if(Month != 12) v1 = 0 ;
	   	    else v1 = K ;
	   	    if(Date != 31) v2 = 0 ;
	   	    else v2 = K;
	   	    f[1][1] = Math.min(v1 , v2) ;

	   	    if(Month == 2) v1 = 0 ;
	   	    else v1 = Math.abs(Month-2)+K ;
	   	    if(Date == 28) v2 = 0 ;
	   	    else v2 = Math.abs(Date-28)+K;
	   	    f[2][0] = v1 + v2 ;
	   	     
	   	    if(Month != 2) v1 = 0 ;
	   	    else v1 = K ;
	   	    if(Date != 28) v2 = 0 ;
	   	    else v2 = K;
	   	    f[2][1] = Math.min(v1 , v2) ;
	   	     
	   	    if(Year%4==0) v1 = 0;
		    else v1 = Math.abs(Year%4-0)+K;
	   	    if(Year%100!=0)v2 = 0 ;
	   	    else v2 = K ;
	   	    if(Year%400==0)v3 = 0 ;
	   	    else v3 = Math.abs(Year%400)+K ;
		    f[3][0] = Math.min(v1+v2, v3) ;
		     
		    if(Year%4!=0) v1 = 0;
		    else v1 = K;
	   	    if(Year%100==0)v2 = 0 ;
	   	    else v2 = Math.abs(Year%100)+K ;
	   	    if(Year%400!=0)v3 = 0 ;
	   	    else v3 = K ;
		    f[3][1] = Math.min(v1,v2)+v3 ;
	   	    
		    if(Month != 12) v1 = 0 ;
	   	    else v1 = K ;
	   	    if(Date == 31) v2 = 0 ;
	   	    else v2 = Math.abs(Date-31)+K ; 
	   	    if(Month == 2) v3 = 0 ;
	   	    else v3 = Math.abs(Month-2)+K ;
	   	    if(Date == 29) v4 = 0 ;
	   	    else v4 = Math.abs(Date-29)+K;
	   	    if(Month == 4) v5 = 0 ;
	   	    else v5 = Math.abs(Month-4)+K ;
	   	    if(Month == 6) v6 = 0 ;
	   	    else v6 = Math.abs(Month-6)+K ;
	   	    if(Month == 9) v7 = 0 ;
	   	    else v7 = Math.abs(Month-9)+K ;
	   	    if(Month == 11) v8 = 0 ;
	   	    else v8 = Math.abs(Month-11)+K ;
	   	    if(Date == 30) v9 = 0 ;
	   	    else v9 = Math.abs(Date-30)+K;
	   	    v10 = v1 + v2 ;
	   	    v11 = v3 + v4 ;
	   	    v12 = Math.min(v5, v6);
	   	    v13 = Math.min(v7, v8);
	   	    v14 = Math.min(v12, v13) + v9;
	   	    v15 = Math.min(v10, v11);
	   	    f[4][0] = Math.min(v15, v14) ;
	   	    
	   	    if(Month == 12) v1 = 0 ;
	   	    else v1 = Math.abs(Month-12)+K ;
	   	    if(Date != 31) v2 = 0 ;
	   	    else v2 = K ; 
	   	    if(Month != 2) v3 = 0 ;
	   	    else v3 = K ;
	   	    if(Date != 29) v4 = 0 ;
	   	    else v4 = K;
	   	    if(Month != 4) v5 = 0 ;
	   	    else v5 = K ;
	   	    if(Month != 6) v6 = 0 ;
	   	    else v6 = K ;
	   	    if(Month != 9) v7 = 0 ;
	   	    else v7 = K ;
	   	    if(Month != 11) v8 = 0 ;
	   	    else v8 = K ;
	   	    if(Date != 30) v9 = 0 ;
	   	    else v9 = K;
	   	    v10 = Math.min(v1, v2);
	   	    v11 = Math.min(v3, v4);
	   	    v12 = Math.min(v5+v6+v7+v8 , v9);
	   	    f[4][1] = v10 + v11 + v12 ;
      	}
	 	if(func_num == 25)  //calculator
	 	{
	 		int ch2 = x[1] ;
      		char cmd = (char) (ch2) ;
      		
      		if(cmd == '+') f[0][0] = 0 ;
      		else f[0][0] = Math.abs(cmd - '+') + K;
      		if(cmd == '-') f[0][1] = 0 ;
      		else f[0][1] = Math.abs(cmd - '-') + K;
      		if(cmd == '.') f[0][2] = 0 ;
      		else f[0][2] = Math.abs(cmd - '.') + K;
      		if(cmd == '/') f[0][3] = 0 ;
      		else f[0][3] = Math.abs(cmd - '/') + K;
      		if(cmd == 'p') f[0][4] = 0 ;
      		else f[0][4] = Math.abs(cmd - 'p') + K;
      		if(cmd == 'a') f[0][5] = 0 ;
      		else f[0][5] = Math.abs(cmd - 'a') + K;
      		if(cmd == 'b') f[0][6] = 0 ;
      		else f[0][6] = Math.abs(cmd - 'b') + K;
      		if(cmd == 'c') f[0][7] = 0 ;
      		else f[0][7] = Math.abs(cmd - 'c') + K;
      		if(cmd == 'd') f[0][8] = 0 ;
      		else f[0][8] = Math.abs(cmd - 'd') + K;
      		if(cmd == 'e') f[0][9] = 0 ;
      		else f[0][9] = Math.abs(cmd - 'e') + K;
      		if(cmd == 'f') f[0][10] = 0 ;
      		else f[0][10] = Math.abs(cmd - 'f') + K;
      		if(cmd == 'S') f[0][11] = 0 ;
      		else f[0][11] = Math.abs(cmd - 'S') + K;
      		if(cmd == 'C') f[0][12] = 0 ;
      		else f[0][12] = Math.abs(cmd - 'C') + K;
      		if(cmd == 'T') f[0][13] = 0 ;
      		else f[0][13] = Math.abs(cmd - 'T') + K;
      		if(cmd!='+' && cmd!='-' && cmd!='.' && cmd!='/' && cmd!='p' && cmd!='a' && cmd!='b' && cmd!='c' && cmd!='d' && cmd!='e' && cmd!='f' && cmd!='S'&& cmd!='C' && cmd!='T' ) 
      			f[0][14] = 0;
      		else
      			f[0][14] = K;
	 	}
	 	if(func_num == 26) //commission
		{
			int totallocks = x[0] ;
			int totalstocks = x[1] ;
			int totalbarrels = x[2] ;
			
			double  lockprice = 45.0 ;
			double  stockprice = 30.0 ;
			double  barrelprice = 25.0 ;
			
			double  locksales = lockprice * totallocks ;
			double  stocksales = stockprice * totalstocks ;
			double  barrelsales = barrelprice * totalbarrels ;
			double  sales = locksales + stocksales + barrelsales ;
			double v1,v2 ;

			if(sales > 1800.0) {f[0][0] = 0 ; f[0][1] = (sales - 1800.0) + K ;}
			else {f[0][0] = (1800.0-sales)+K ; f[0][1] = 0 ;}
			 
			if( sales >500.0)v1 = 0 ;
			else v1 = (500.0-sales) + K ;
			if(sales <= 1800.0)v2 = 0 ;
			else v2 = (sales-1800.0) +K ;
			f[1][0] = v1 + v2 ;
			
			if( sales <=500.0)v1 = 0 ;
			else v1 = (sales-500.0) + K ;
			if(sales > 1800.0)v2 = 0 ;
			else v2 = (1800.0-sales) + K ;
			f[1][1] = Math.min(v1, v2) ;		
		}
	 	if(func_num == 27) //premium
		{
			int  driverage = x[0] ;
			int  points = x[1] ;
			double v1,v2,v3,v4,v5,v6,v7,v8,v9,v10 ;	
			
			if(driverage >=16)v1 = 0 ;
			else  v1 = (16-driverage)+K ;
			if(driverage < 20)v2 = 0 ;
			else v2 = (driverage-20)+K ;
			f[0][0] = v1 + v2 ;
			
			if(driverage < 16)v1 = 0 ;
			else  v1 = (driverage-16)+K ;
			if(driverage >= 20)v2 = 0 ;
			else v2 = (20-driverage)+K ;
			f[0][1] = Math.min(v1, v2) ;
			
			if(points <= 1){f[1][0] = 0 ; f[1][1] = (1-points)+K ; ;}
			else {f[1][0] = (points-1)+K ; f[1][1] = 0 ;}
			
			if(driverage >= 20)v3 = 0 ;
			else v3 = (20-driverage)+K ;
			if(driverage < 25) v4 = 0 ;
			else v4 = (driverage-25)+K ;
			f[2][0] = v3 + v4 ;
			
			if(driverage < 20)v3 = 0 ;
			else v3 = (driverage-20)+K ;
			if(driverage >= 25) v4 = 0 ;
			else v4 = (25-driverage)+K ;
			f[2][1] = Math.min(v3 , v4) ;
			
			if(points < 3){f[3][0] = 0 ; f[3][1] = (3-points)+K ;}
			else {f[3][0] = (points-3)+K ; f[3][1] = 0 ;}
			
			if(driverage >= 25)v5 = 0 ;
			else v5 = (25-driverage)+K ;
			if(driverage < 45) v6 = 0 ;
			else v6 = (driverage-45)+K ;
			f[4][0] = v5 + v6 ;
			
			if(driverage < 25)v5 = 0 ;
			else v5 = (driverage-25)+K ;
			if(driverage >= 45) v6 = 0 ;
			else v6 = (45-driverage)+K ;
			f[4][1] = Math.min(v5 , v6) ;
			
			if(points < 5){f[5][0] = 0 ; f[5][1] = (5-points)+K ;}
			else {f[5][0] = (points-5)+K ; f[5][1] = 0 ;}
			
			if(driverage >= 45)v7 = 0 ;
			else v7 = (45-driverage)+K ;
			if(driverage < 60) v8 = 0 ;
			else v8 = (driverage-60)+K ;
			f[6][0] = v7 + v8 ;
			
			if(driverage < 45)v7 = 0 ;
			else v7 = (driverage-45)+K ;
			if(driverage >= 60) v8 = 0 ;
			else v8 = (60-driverage)+K ;
			f[6][1] = Math.min(v7 , v8) ;
			
			if(points < 7){f[7][0] = 0 ; f[7][1] = (7-points)+K ;}
			else {f[7][0] = (points-7)+K ; f[7][1] = 0 ;}
			
			if(driverage >= 60)v9 = 0 ;
			else v9 = (60-driverage)+K ;
			if(driverage < 100) v10 = 0 ;
			else v10 = (driverage-100)+K ;
			f[8][0] = v9 + v10 ;
			
			if(driverage < 60)v9 = 0 ;
			else v9 = (driverage-60)+K ;
			if(driverage >= 100) v10 = 0 ;
			else v10 = (100-driverage)+K ;
			f[8][1] = Math.min(v9 , v10) ;
			
			if(points < 5){f[9][0] = 0; f[9][1] = (5-points)+K ;}
			else {f[9][0] = (points-5)+K ; f[9][1] = 0 ;}			
		}
	 	
	 	if(func_num == 28) //decision
		{
			int a = x[0] ;
			int b = x[1] ;
			int c = x[2] ;	
			int d = x[3] ;
					
			if(a==3971)	f[0][0] = 0;
			else	f[0][0] = Math.abs(a-3971)+K;
			if(a==5085)	f[0][1] = 0;
			else	f[0][1] = Math.abs(a-5085)+K;
			if(a==5174)	f[0][2] = 0;
			else	f[0][2] = Math.abs(a-5174)+K;
			if(a!=3971&&a!=5085&&a!=5174)	f[0][3] = 0;
			else	f[0][3] = 3*K;
			
			if(c==5448) f[1][0] = 0;
			else	f[1][0] = Math.abs(c-5448)+K;
			if(c==2463)	f[1][1] = 0;
			else	f[1][1] = Math.abs(c-2463)+K;
			if(c!=5448&&c!=2463)	f[1][2] = 0;
			else	f[1][2] = 2*K;
			
			if(b==4040)	f[2][0] = 0;
			else	f[2][0] = Math.abs(b-4040)+K;
			if(b==5448) f[2][1] = 0;
			else	f[2][1] = Math.abs(b-5448)+K;
			if(b==3268) f[2][2] = 0;
			else	f[2][2] = Math.abs(b-3268)+K;
			if(b!=4040&&b!=5448&&b!=3268)	f[2][3] = 0;
			else	f[2][3] = 3*K;
			
			if(d==5148) f[3][0] = 0;
			else	f[3][0] = Math.abs(d-5148)+K;
			if(d==4662) f[3][1] = 0;
			else	f[3][1] = Math.abs(d-4662)+K;
			if(d!=5148&&d!=4662)	f[3][2] = 0;
			else	f[3][2] = 2*K;
		}

		if(path_num == -1)          //没有目标路径的情况
		{
			for(int k = 0 ; k < NODENUM ; k++)
	   		{
		   		if(visit[k][0] && visit[k][1])
		   			fit[k] = 0 ;
		   		else if(visit[k][0] && (!visit[k][1]))
		   			fit[k] = 1/(f[k][1] + alpha) ;
		   		else if((!visit[k][0]) && visit[k][1])
		   			fit[k] = 1/(f[k][0] + alpha) ;
		   		else
		   			fit[k] = 1/alpha ;
	   		}
		}
		else {// 存在目标路径的情况
			for (int k = 0; k < NODENUM; k++) {
				if (PATH[path_num].charAt(k) == ' ')
					fit[k] = 0;
				else
					fit[k] = 1 / (f[k][PATH[path_num].charAt(k) - '0'] + alpha);

			}

		}
		for(int i=0;i<NODENUM;i++)
			Fitness += fit[i];
		
		return Fitness;
	}
	static boolean isRun(int year) //To identify whether the year is bissextile.  
	 {  
	      if((year%4==0 && year%100!=0) || (year%400==0))  
	      {  
	          return true;  
	      }  
	     else  
	     {  
	        return false;  
	     }  
	 } 
	
	/*返回适应值从高到低排序对应的序号*/
	  public static int[] selectsort (double[] fitness)
	  {
			  int[] sortnum = new int[fitness.length] ;
			  for(int i =0 ; i < fitness.length ; i++)
				  sortnum[i]=i;
			  int max = 0;
			  double tmp = 0;
			  int tmp2 = 0 ;
			  for(int i=0;i<fitness.length;i++)
			  {
			       max = i;
			       /**查找第 i大的数，直到记下第 i大数的位置***/
			       for(int j=i+1;j<fitness.length;j++)
			       {
			            if(fitness[max]<fitness[j]) 
			            max = j;//记下较大数位置，再次比较，直到最大
			       }
			        /***如果第 i大数的位置不在 i,则交换****/
			        if(i!=max)
			        {
					    tmp = fitness[i];
					    fitness[i] = fitness[max];
					    fitness[max] = tmp;
					    
					    tmp2 = sortnum[i];
					    sortnum[i] = sortnum[max];
					    sortnum[max] = tmp2;
			        }
			  }
			  return sortnum ;
	  }
	 
	//获取平均值
	  static double getAverage(int[] array , int num){
	      int sum = 0;
	      for(int i = 0;i < num;i++){
	          sum += array[i];
	      }
	      return (double)(sum / num);
	  }
	 
	  //标准差
	  static double getStandardDevition(int[] array , int num){
	      double sum = 0;
	      for(int i = 0;i < num;i++){
	          sum += Math.sqrt(((double)array[i] -getAverage(array, num)) * (array[i] -getAverage(array, num)));
	      }
	      return (sum / (num - 1));
	  } 
	  //返回最大值下标
	  static int getBestIndex(double[]array)
	  {
		  int index;
		  index = 0;
		  for(int i=1;i<array.length;i++)
			  if(array[i]>array[index])
				  index = i;
		  return index;
	  }
//	  static int[] getIndex1(int lb, int ub, int best, double step) {
//			int[] index = new int[step_length1+1];
//			int temp = -1;
//			if (best + step_length1/2 * step <= ub && best - step_length1/2 * step >= lb) {
//				for (int i = 0; i < step_length1+1; i++)
//					index[i] = (int) (step * (i - step_length1/2) + best);
//			} else {
//				for (int i = 1; i < step_length1+1; i++)
//					if (best - i * step < lb) {
//						temp = i;
//						break;
//					}
//				if (temp == -1) {
//					for (int i = 0; i < step_length1+1; i++)
//						index[i] = (int) (best - step * i) ;
//				} else {
//					for (int i = 0; i < step_length1+1; i++)
//						index[i] = (int) (best + step * (-temp + i + 1));
//				}
//			}
//
//			return index;
//		}
	  static int[] getIndex1(int lb, int ub, int best, double step) {
			int[] index = new int[step_length1+1];
			int temp1 = -1, temp2 = -1;
			if (best + step_length1/2 * step <= ub && best - step_length1/2 * step >= lb) {
				for (int i = 0; i < step_length1+1; i++)
					index[i] = (int) (step * (i - step_length1/2) + best);
			} else {
				for (int i = 1; i < step_length1/2+1; i++)
				{
					if (best - i * step < lb) {
						temp1 = i;
						break;
					}
					if(best + i * step > ub){
						temp2 = i;
						break;
					}
				}
				if (temp1 > -1) {
					for (int i = 0; i < step_length1+1; i++)
						index[i] = (int) (best + step * (-temp1 + i + 1));
				} else if(temp2 > -1){
					for (int i = 0; i < step_length1+1; i++)
						index[i] = (int) (best + step * ( temp2 - i - 1));
				}
				else{
					for (int i = 0; i < step_length1+1; i++)
						index[i] = (int) (best - step * i) ;
				}
			}

			return index;
		}
	  static int[] getIndex2(int lb, int ub, int best, double step) {
			int[] index = new int[step_length2+1];
			int temp1 = -1, temp2 = -1;
			if (best + step_length2/2 * step <= ub && best - step_length2/2 * step >= lb) {
				for (int i = 0; i < step_length2+1; i++)
					index[i] = (int) (step * (i - step_length2/2) + best);
			} else {
				for (int i = 1; i < step_length2/2+1; i++)
				{
					if (best - i * step < lb) {
						temp1 = i;
						break;
					}
					if(best + i * step > ub){
						temp2 = i;
						break;
					}
				}
				if (temp1 > -1) {
					for (int i = 0; i < step_length2+1; i++)
						index[i] = (int) (best + step * (-temp1 + i + 1));
				} else if(temp2 > -1){
					for (int i = 0; i < step_length2+1; i++)
						index[i] = (int) (best + step * ( temp2 - i - 1));
				}
				else{
					for (int i = 0; i < step_length2+1; i++)
						index[i] = (int) (best - step * i) ;
				}
			}

			return index;
		}
	  static int[] getIndex3(final int lb, final int ub, final int best, final double step)
	  {
		  final int[] index = new int[step_length1];
		  index[0] = (int) (best - step);
		  index[1] = (int) (best + step);
		  return index;
	  }
	  
	  static int[] createTestCase(int[] input, int dim, int value)
	  {
		  int[] newTestCase = new int[R];
		  for(int i=0;i<input.length;i++){
			  if(i!=dim)
				  newTestCase[i] = input[i];
			  else
				  newTestCase[i] = value;
		  }
		  return newTestCase;
	  }
}
