/*
 * Author: tdanford
 * Date: Dec 19, 2008
 */
package org.seqcode.utils.models;

import java.io.*;

import org.seqcode.utils.Closeable;
import org.seqcode.utils.json.*;


public interface ModelInput<M extends Model> extends Closeable {

	public M readModel();
	
	public static class LineReader<T extends Model> implements ModelInput<T> {
		
		private Class<T> modelClass;
		private BufferedReader br;
		
		public LineReader(Class<T> cls, Reader r) { 
			modelClass = cls;
			br = new BufferedReader(r);
		}

		public LineReader(Class<T> cls, InputStream is) { 
			this(cls, new InputStreamReader(is));
		}

		public T readModel() {
			try {
				String line = br.readLine();
				if(line != null) { 
					JSONObject jsonObject = new JSONObject(new JSONTokener(line));
					Object unjsoned = Model.unjsonify(modelClass, jsonObject);
					
					if(unjsoned != null && 
					   Model.isSubclass(unjsoned.getClass(), modelClass)) {
						
						return (T)unjsoned;
					}
				} else { 
					close();
				}
			} catch (JSONException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}

			return null;
		}

		public void close() {
			if(isClosed()) { return; }
			try {
				br.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			br = null;
		}

		public boolean isClosed() {
			return br == null;
		} 
	}
}
